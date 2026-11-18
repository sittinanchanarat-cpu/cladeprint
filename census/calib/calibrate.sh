#!/usr/bin/env bash
# Coverage calibration: how scores and calls degrade with sequencing depth.
#
# Per run (census/calib_runs.tsv: well-covered literature-labelled training runs, one per genotype
# cluster, 4 per species / lineage, whose leave-one-out sensitivity is >= 0.95 so being in the
# training set barely helps them): fetch from the AWS SRA copy, decode a 4.2M-spot hash subsample,
# then nested subsamples at DEPTHS spots - every depth uses the same seed at an absolute fraction
# (depth / run spots), so each smaller sample is a subset of the larger ones. The seed differs from
# the census subsample (20261012), so the reads are not the ones the bundle was trained on.
# Each depth is classified against BUNDLE (KMC counting inside classify); reads, databases and the
# .sra are deleted after use. Results: OUT/<run>.tsv (classify --long, one row per depth, sample
# name <run>@<spots>), OUT/<run>.spots (depth, spots kept).
#
# Usage: micromamba run -n cladeprint bash census/calib/calibrate.sh
#   JOBS=2 runs at once; LIMIT=N first N runs; RUNS="SRRx SRRy" only these; OUT=census/calib/out;
#   BUNDLE=census/out/bundle; DECODE_JOBS=4; THREADS=4 per classify; HEADROOM_MB=6144
set -euo pipefail
ROOT=$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)
cd "$ROOT"
OUT=${OUT:-census/calib/out}
BUNDLE=${BUNDLE:-census/out/bundle}
SEED=20261020
TOP=4200000
DEPTHS="125000 250000 500000 1000000 2000000 4000000"
JOBS=${JOBS:-2}
THREADS=${THREADS:-4}
HEADROOM_MB=${HEADROOM_MB:-6144}  # free space beyond the .sra before a run starts
LIMIT=${LIMIT:-1000}
mkdir -p $OUT/tmp
source scripts/sra.sh  # fetch, decode, CLADEPRINT
for t in fastq-dump kmc kmc_tools curl "$CLADEPRINT"; do
  command -v "$t" > /dev/null || { echo "FAIL: $t not found (run inside the cladeprint env)"; exit 1; }
done

one() {
  local run=$1 spots=$2 size_mb=$3
  [[ -s $OUT/$run.tsv ]] && return 0
  local T=$OUT/tmp/$run
  mkdir -p $T
  # wait for disk: the .sra, ~3 GB of reads at all depths, and KMC temporary space
  while (( $(df -Pm $OUT | awk 'NR==2{print $4}') < ${size_mb%.*} + HEADROOM_MB )); do sleep 60; done
  local url="https://sra-pub-run-odp.s3.amazonaws.com/sra/$run/$run"
  wait_net "$url"
  local size; size=$(curl -sSI --retry 5 "$url" | tr -d '\r' | awk 'tolower($1)=="content-length:"{print $2}')
  [[ $size =~ ^[0-9]+$ ]] || { echo "FAIL $run (not in the AWS copy)"; return 1; }
  local tries=0
  until [[ $(stat -c %s $T/$run.sra 2>/dev/null || echo 0) -eq $size ]]; do
    (( ++tries <= 20 )) || { echo "FAIL $run (download incomplete)"; return 1; }
    wait_net "$url"
    fetch "$url" $size $T/$run.sra || true  # short parts are resumed on the next try
  done
  local frac seen kept
  frac=$(awk -v n=$TOP -v s=$spots 'BEGIN{f=n/s; print (f>1?1:f)}')
  read seen kept < <(decode $T/$run.sra $spots $frac $SEED $T/top.fq.gz)
  rm -f $T/$run.sra
  [[ $seen -eq $spots ]] || { echo "FAIL $run (stream had $seen of $spots spots)"; return 1; }
  printf 'sample\treads\n' > $T/sheet.tsv
  : > $T/spots
  local d f
  for d in $DEPTHS; do
    f=$(awk -v n=$d -v s=$spots 'BEGIN{print n/s}')
    zcat $T/top.fq.gz | "$CLADEPRINT" subsample --fraction $f --seed $SEED --counts $T/$d.count \
      | gzip -1 > $T/$d.fq.gz
    echo -e "$d\t$(cut -d' ' -f2 $T/$d.count)" >> $T/spots
    printf '%s@%s\t%s\n' $run $d $d.fq.gz >> $T/sheet.tsv  # relative to the sheet's directory
  done
  rm -f $T/top.fq.gz
  "$CLADEPRINT" --threads $THREADS classify $BUNDLE --samples $T/sheet.tsv --long \
    --tmp $T/kmc --memory-gb 4 -o $T/$run.tsv 2> $T/classify.log \
    || { echo "FAIL $run (classify; log kept in $T/classify.log)"; return 1; }
  local rows; rows=$(tail -n +2 $T/$run.tsv | wc -l)
  [[ $rows -eq $(wc -w <<< "$DEPTHS") ]] \
    || { echo "FAIL $run (classify wrote $rows rows; log kept in $T/classify.log)"; return 1; }
  mv $T/spots $OUT/$run.spots
  mv $T/$run.tsv $OUT/$run.tsv  # written last: marks the run complete
  rm -rf $T
  echo "done $run"
}
# each worker is a fresh `bash -c` that does not inherit set -e, and errexit is ignored inside any
# `cmd || ...`, so `one` runs in its own errexit subshell and the exit code is checked afterwards
export -f one; export OUT BUNDLE SEED TOP DEPTHS THREADS HEADROOM_MB

if [[ -n ${RUNS:-} ]]; then
  tail -n +2 census/calib_runs.tsv | awk -F'\t' -v r=" $RUNS " 'index(r, " "$1" ") {print $1, $7, $8}'
else
  tail -n +2 census/calib_runs.tsv | awk -F'\t' '{print $1, $7, $8}' | head -n $LIMIT
fi | xargs -P $JOBS -L1 bash -c \
  '( set -euo pipefail; one "$@" ); rc=$?; (( rc == 0 )) || echo "FAIL $1 (exit $rc)"' _ | tee -a $OUT/calibrate.log
echo "finished: $(ls $OUT/*.tsv 2>/dev/null | wc -l) runs"
