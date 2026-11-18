#!/usr/bin/env bash
# Coverage-matched recount of the census runs whose 1M-spot database has no usable genome peak
# (low_coverage = true in out/db/RUN.stats.toml: mostly old, error-rich uploads, where many reads
# fail the species screen's quality filter). Each run is subsampled again with the same seed; the
# sampler keeps a spot by a hash of its name, so a larger sample contains the 1M-spot one. The spot
# count aims at the census median peak (TARGET=6): n = n * TARGET / peak when a peak was found,
# else n * 3, capped at all spots; up to MAX_TRIES attempts on the same download. A run is accepted
# when its peak is in [LO, HI] (5-12, the bulk of the census) or all spots are used; runs still
# without a peak stay low_coverage and are left out by cladeprint's QC.
#
# The 1M-spot outputs move to out/db_1M/ (RUN.kmc_*, .hist, .stats.toml, .species.tsv) for the
# before/after table. One line per attempt in out/logs/recount.tsv. Resumes: a run is skipped once
# out/db/RUN.recount_done exists.
#
# Usage (after fetch_count.sh has finished): micromamba run -n cladeprint bash census/recount.sh
#   JOBS=4 KMC_THREADS=4 as fetch_count.sh; TARGET=6 LO=5 HI=12 MAX_TRIES=4
set -euo pipefail
cd "$(dirname "$0")"
source fetch_count.sh  # one, OUT, NSPOTS, ...
TARGET=${TARGET:-6} LO=${LO:-5} HI=${HI:-12} MAX_TRIES=${MAX_TRIES:-4}
mkdir -p $OUT/db_1M
LOG=$OUT/logs/recount.tsv
[[ -s $LOG ]] || echo -e "run\ttry\tnspots\tspots_total\tkept\tcutoff\tpeak\tsolid_kmers\terror_kmers\tlow_coverage" > $LOG

stat_of() { awk -v k=$2 '$1 == k {print $3}' $1.stats.toml; }  # db key

rec() {
  run=$1 set=$2 spots=$3 size_mb=$4
  local db=$OUT/db/$run old=$OUT/db_1M/$run f n p try out
  [[ -e $db.recount_done ]] && return 0
  if [[ ! -s $old.stats.toml ]]; then
    for f in $db.kmc_pre $db.kmc_suf $db.hist; do [[ -e $f ]] && mv $f $OUT/db_1M/; done
    cp $OUT/species/$run.tsv $old.species.tsv
    mv $db.stats.toml $old.stats.toml  # last: marks the move as complete
  fi
  p=$(stat_of $old peak)
  n=$(( p > 0 ? NSPOTS * TARGET / p : NSPOTS * 3 ))
  for ((try = 1; try <= MAX_TRIES; try++)); do
    (( n > spots )) && n=$spots
    rm -f $db.kmc_pre $db.kmc_suf $db.hist $db.stats.toml
    out=$(NSPOTS=$n KEEP_SRA=1 one $run $set $spots $size_mb) || { echo "FAIL $run (recount try $try): $out"; return 1; }
    p=$(stat_of $db peak)
    echo -e "$run\t$try\t$n\t$spots\t$(awk '{print $3}' <<< "$out")\t$(stat_of $db cutoff)\t$p\t$(stat_of $db solid_kmers)\t$(stat_of $db error_kmers)\t$(stat_of $db low_coverage)" >> $LOG
    (( (p >= LO && p <= HI) || n >= spots )) && break
    local next=$(( p > 0 ? n * TARGET / p : n * 3 ))
    (( next == n )) && break
    n=$next
  done
  rm -f $OUT/tmp/$run.sra $OUT/tmp/${run}{,_1,_2}.fastq.gz
  touch $db.recount_done
  echo "recounted $run: try $try, $n of $spots spots, peak $p, low_coverage $(stat_of $db low_coverage)"
}
export -f rec stat_of; export TARGET LO HI MAX_TRIES LOG

# runs to recount: low coverage at 1M spots (still in db/, or already moved to db_1M/)
shopt -s nullglob
low=$(for f in $OUT/db/*.stats.toml $OUT/db_1M/*.stats.toml; do
  if grep -q '^low_coverage = true' $f; then basename $f .stats.toml; fi; done | sort -u)
echo "$(wc -w <<< "$low") runs to recount"
tail -n +2 ${SHEET:-census_runs.tsv} \
  | awk -F'\t' 'NR == FNR {l[$1] = 1; next} ($2 in l) {print $2, $3, $11, $12}' <(tr ' ' '\n' <<< "$low") - \
  | xargs -P $JOBS -L1 bash -c 'rec "$@" || echo "FAIL $1"' _ | tee -a $OUT/logs/recount.log
done_runs=($OUT/db/*.recount_done) still=0
for f in "${done_runs[@]}"; do
  if grep -q '^low_coverage = true' ${f%.recount_done}.stats.toml; then still=$((still + 1)); fi
done
echo "recount finished: ${#done_runs[@]} runs; still low coverage: $still"
