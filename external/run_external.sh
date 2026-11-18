#!/usr/bin/env bash
# External test set: fetch each run -> random 1M-pair subsample -> KMC presence database, exactly as
# for the training isolates (pilot/run_pilot.sh subsample + kmers/build_dbs.sh cutoff), then score.
#
# Downloads use 8 parallel byte ranges per run (S3 throttles single connections to ~0.4 MB/s here).
# Differences from training, both for disk: reads are streamed out of the .sra (fastq-dump --stdout)
# straight into the subsample, so no full FASTQ is written (fasterq-dump is not used: even with -Z it
# writes the whole FASTQ to its temp directory, ~5x the .sra, which is what failed the large runs);
# no H99 mapping (training BAMs held every read, mapped or not, so the k-mer input is the same). The stream must contain every spot in the
# run (runinfo count), which catches truncated downloads. A run is skipped once db/RUN exists.
#
# Usage: micromamba run -n cladeprint bash run_external.sh      (JOBS=2, LIMIT=N for the first N runs,
#        SET=pure|hybrid|reserve to take one set, SHEET=file for another sheet; overnight.sh drives it)
set -euo pipefail
cd "$(dirname "$0")"
K=31
NPAIRS=1000000
SEED=20261012  # same seed as pilot/run_pilot.sh
JOBS=${JOBS:-2}
LIMIT=${LIMIT:-100000}
mkdir -p tmp db hist logs
[[ -f stats.tsv ]] || echo -e "run\tspots\tspots_seen\tpairs_kept\tcutoff\tpeak_cov\tsolid_kmers\terror_kmers\tlow_cov" > stats.tsv

source ../scripts/sra.sh  # fetch, decode
[[ -x $CLADEPRINT ]] || { echo "build cladeprint first (cargo build --release)" >&2; exit 1; }

one() {
  run=$1 spots=$2 size_mb=$3
  [[ -s db/$run.kmc_suf ]] && return 0
  # wait for disk: this run's .sra plus 4 GB of headroom (two jobs may hold large .sra files at once)
  while (( $(df -Pm . | awk 'NR==2{print $4}') < ${size_mb%.*} + 4096 )); do sleep 60; done
  t0=$SECONDS
  url="https://sra-pub-run-odp.s3.amazonaws.com/sra/$run/$run"
  size=$(curl -sSI --retry 5 "$url" | tr -d '\r' | awk 'tolower($1)=="content-length:"{print $2}')
  [[ $size =~ ^[0-9]+$ ]] || { echo "FAIL $run (not on AWS)"; return 1; }
  [[ $(stat -c %s tmp/$run.sra 2>/dev/null || echo 0) -eq $size ]] || fetch "$url" $size tmp/$run.sra  # reuse a complete download
  [[ $(stat -c %s tmp/$run.sra 2>/dev/null || echo 0) -eq $size ]] || { rm -f tmp/$run.sra; echo "FAIL $run (download incomplete)"; return 1; }

  frac=$(awk -v n=$NPAIRS -v s=$spots 'BEGIN{f=1.05*n/s; print (f>1?1:f)}')
  t_dl=$SECONDS
  read seen kept < <(decode tmp/$run.sra $spots $frac $SEED tmp/$run.fq.gz)
  rm -f tmp/$run.sra
  [[ $seen -eq $spots ]] || { rm -f tmp/$run.fq.gz; echo "FAIL $run (stream had $seen of $spots spots)"; return 1; }

  t_kmc=$SECONDS
  # presence database: same steps and cutoff rule as kmers/build_dbs.sh
  mkdir -p tmp/$run.kmc
  kmc -k$K -m2 -t3 -ci1 -cs1000 tmp/$run.fq.gz tmp/$run.all tmp/$run.kmc/ > /dev/null 2>&1
  rm -rf tmp/$run.fq.gz tmp/$run.kmc
  kmc_tools transform tmp/$run.all histogram hist/$run.txt -cx1000 > /dev/null 2>&1
  read cut peak solid err low < <(awk '{h[$1]=$2} END{
      c=2; while (c<20 && h[c+1]<h[c]) c++;
      low=0; if (c>=20) {c=2; low=1}
      p=c; for (i=c; i<=1000; i++) if (h[i]>h[p]) p=i;
      if (low) p=0; else if (p<5) low=1;
      s=0; for (i=c; i<=1000; i++) s+=h[i]; e=0; for (i=1; i<c; i++) e+=h[i];
      print c, p, s, e, low}' hist/$run.txt)
  kmc_tools transform tmp/$run.all -ci$cut set_counts 1 tmp/$run.db > /dev/null 2>&1
  rm -f tmp/$run.all.kmc_pre tmp/$run.all.kmc_suf
  echo "$run download $((t_dl - t0)) s, stream+subsample $((t_kmc - t_dl)) s, kmc $((SECONDS - t_kmc)) s" >> logs/timing.log
  mv tmp/$run.db.kmc_pre db/$run.kmc_pre && mv tmp/$run.db.kmc_suf db/$run.kmc_suf  # db appears only when complete
  echo -e "$run\t$spots\t$seen\t$kept\t$cut\t$peak\t$solid\t$err\t$low" >> stats.tsv
  echo "done $run: $kept of $seen spots, cutoff $cut, peak ${peak}x, $solid solid k-mers"
}
export -f one; export K NPAIRS SEED

SHEET=${SHEET:-samples.tsv}  # any sheet with samples.tsv columns; SET=pure|hybrid|reserve filters it
tail -n +2 $SHEET | awk -F'\t' -v set="${SET:-}" 'set == "" || $8 == set {print $5, $6, $7}' | head -n $LIMIT \
  | xargs -P $JOBS -L1 bash -c 'one "$@" || echo "FAIL $1"' _ | tee -a logs/fetch.log
echo "fetch finished: $(ls db/*.kmc_suf 2>/dev/null | wc -l) databases"
