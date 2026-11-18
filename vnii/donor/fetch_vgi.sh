#!/usr/bin/env bash
# Donor search: fetch each VGI run in vgi_samples.tsv (or SHEET=file) -> random 1M-spot subsample -> map to H99 + WM276
# (vnii/locus/ref2.mmi) -> keep only the regions needed, so each run leaves a few MB on disk:
#   the WM276 block window and the H99 locus window, and a 1 Mb baseline on each genome
#   (WM276 chr A, H99 chr 1) for genome-wide distance and the species check.
# Download and subsampling as in external/run_external.sh (8 byte ranges from the AWS SRA copy, whole
# spots kept with the same seed, stream must contain every spot). Single-end runs map as single reads;
# for paired runs minimap2 pairs consecutive reads with the same name. A run is skipped once bam/RUN.bam exists.
#
# Usage: micromamba run -n cladeprint bash vnii/donor/fetch_vgi.sh     (JOBS=2)
set -euo pipefail
cd "$(dirname "$0")"
NSPOTS=1000000
SEED=20261012  # same seed as pilot/run_pilot.sh
JOBS=${JOBS:-2}
REGIONS="gattii|NC_014943.1:900000-970000 neoformans|NC_026751.1:960000-1050000 gattii|NC_014938.1:200000-1200000 neoformans|NC_026745.1:200000-1200000"
mkdir -p tmp bam logs

source ../../scripts/sra.sh  # fetch, decode
[[ -x $CLADEPRINT ]] || { echo "build cladeprint first (cargo build --release)" >&2; exit 1; }

one() {
  run=$1 spots=$2 size_mb=$3
  [[ -s bam/$run.bam ]] && return 0
  while (( $(df -Pm . | awk 'NR==2{print $4}') < ${size_mb%.*} + 4096 )); do sleep 60; done
  url="https://sra-pub-run-odp.s3.amazonaws.com/sra/$run/$run"
  size=$(curl -sSI --retry 5 "$url" | tr -d '\r' | awk 'tolower($1)=="content-length:"{print $2}')
  [[ $size =~ ^[0-9]+$ ]] || { echo "FAIL $run (not on AWS)"; return 1; }
  [[ $(stat -c %s tmp/$run.sra 2>/dev/null || echo 0) -eq $size ]] || fetch "$url" $size tmp/$run.sra  # reuse a complete download
  [[ $(stat -c %s tmp/$run.sra 2>/dev/null || echo 0) -eq $size ]] || { rm -f tmp/$run.sra; echo "FAIL $run (download incomplete)"; return 1; }

  frac=$(awk -v n=$NSPOTS -v s=$spots 'BEGIN{f=1.05*n/s; print (f>1?1:f)}')
  read seen kept < <(decode tmp/$run.sra $spots $frac $SEED tmp/$run.fq.gz)
  rm -f tmp/$run.sra
  [[ $seen -eq $spots ]] || { rm -f tmp/$run.fq.gz; echo "FAIL $run (stream had $seen of $spots spots)"; return 1; }

  minimap2 -t2 -ax sr ../locus/ref2.mmi tmp/$run.fq.gz 2>/dev/null | samtools sort -@1 -m 400M -o tmp/$run.all.bam -
  samtools index tmp/$run.all.bam
  samtools view -b -o tmp/$run.bam tmp/$run.all.bam $REGIONS
  echo -e "$run\t$(samtools view -c -F 0x904 tmp/$run.all.bam)\t$(samtools view -c -F 0x904 tmp/$run.all.bam neoformans\|NC_026745.1 neoformans\|NC_026746.1)\t$(samtools view -c -F 0x904 tmp/$run.all.bam gattii\|NC_014938.1 gattii\|NC_014939.1)" >> logs/mapped.tsv
  rm -f tmp/$run.fq.gz tmp/$run.all.bam tmp/$run.all.bam.bai
  mv tmp/$run.bam bam/$run.bam && samtools index bam/$run.bam  # bam appears only when complete
  echo "done $run: $kept of $seen spots"
}
export -f one; export NSPOTS SEED REGIONS

tail -n +2 ${SHEET:-vgi_samples.tsv} | awk -F'\t' '{print $4, $5, $6}' \
  | xargs -P $JOBS -L1 bash -c 'one "$@" || echo "FAIL $1"' _ | tee -a logs/fetch.log
echo "fetch finished: $(ls bam/*.bam 2>/dev/null | wc -l) runs"
