#!/usr/bin/env bash
# Census, per run (census/census_runs.tsv): fetch from the AWS SRA copy -> random 1M-spot subsample
# (parallel decode, scripts/sra.sh) -> species screen on ~50k spots spread over the subsample
# (competitive mapping to the panel of the run's set: the 7 species of the Cryptococcus complex, as
# pilot/species_screen.sh, or 6 Sporothrix species) -> presence
# database with `cladeprint count` (KMC, histogram cutoff; writes <db>.stats.toml last). The reads
# are deleted after use. A run is skipped once out/db/RUN.stats.toml exists, so the script resumes.
#
# Outputs under OUT (default census/out): db/RUN.kmc_pre|kmc_suf|hist|stats.toml, species/RUN.tsv,
# logs/. Species columns follow PANEL_DIR/<panel>.species.txt. Needs: curl, fastq-dump, minimap2,
# samtools, gzip, kmc, kmc_tools, cladeprint (built), the species panels (pilot/fetch_refs.sh:
# PANEL_DIR/{complex,sporothrix}.{mmi,species.txt,nuclear.txt}).
#
# Usage: micromamba run -n cladeprint bash census/fetch_count.sh
#   JOBS=6 runs at once; DECODE_JOBS=4 fastq-dump per run; KMC_THREADS=4; SHEET=census_runs.tsv;
#   SET=cryptococcus|sporothrix (default both); LIMIT=N first N runs; OUT=out; NSPOTS=1000000 spots
#   per subsample; KEEP_SRA=1 keeps the downloaded .sra (census/recount.sh). Runs missing from the AWS
#   copy are fetched as FASTQ from ENA, mates interleaved, and subsampled the same way.
#   Sourced (census/recount.sh), it only defines `one` and does not start the census.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
OUT=${OUT:-out}
NSPOTS=${NSPOTS:-1000000}
SEED=20261012
SCREEN_SPOTS=50000
JOBS=${JOBS:-6}
KMC_THREADS=${KMC_THREADS:-4}
LIMIT=${LIMIT:-100000}
PANEL_DIR=${PANEL_DIR:-../pilot/ref/panel}
mkdir -p $OUT/db $OUT/species $OUT/logs $OUT/tmp
source ../scripts/sra.sh  # fetch, decode, CLADEPRINT
[[ -x $CLADEPRINT ]] || { echo "build cladeprint first (cargo build --release)" >&2; exit 1; }
for p in complex sporothrix; do
  for f in mmi species.txt nuclear.txt; do
    [[ -s $PANEL_DIR/$p.$f ]] || { echo "species panel $PANEL_DIR/$p.$f missing (pilot/fetch_refs.sh)" >&2; exit 1; }
  done
done

one() {
  run=$1 set=$2 spots=$3 size_mb=$4
  local panel=$PANEL_DIR/$([[ $set == sporothrix ]] && echo sporothrix || echo complex)
  local T=$OUT/tmp db=$OUT/db/$run
  [[ -s $db.stats.toml ]] && return 0
  # wait for disk: this run's .sra plus 4 GB of headroom
  while (( $(df -Pm $OUT | awk 'NR==2{print $4}') < ${size_mb%.*} + 4096 )); do sleep 60; done
  local t0=$SECONDS
  local url="https://sra-pub-run-odp.s3.amazonaws.com/sra/$run/$run"
  local size; size=$(curl -sSI --retry 5 "$url" | tr -d '\r' | awk 'tolower($1)=="content-length:"{print $2}')
  local frac; frac=$(awk -v n=$NSPOTS -v s=$spots 'BEGIN{f=1.05*n/s; print (f>1?1:f)}')
  local t_dl seen kept
  if [[ $size =~ ^[0-9]+$ ]]; then
    [[ $(stat -c %s $T/$run.sra 2>/dev/null || echo 0) -eq $size ]] || fetch "$url" $size $T/$run.sra
    [[ $(stat -c %s $T/$run.sra 2>/dev/null || echo 0) -eq $size ]] || { rm -f $T/$run.sra; echo "FAIL $run (download incomplete)"; return 1; }
    t_dl=$SECONDS
    read seen kept < <(decode $T/$run.sra $spots $frac $SEED $T/$run.fq.gz)
    [[ ${KEEP_SRA:-0} == 1 ]] || rm -f $T/$run.sra
  else
    # not in the AWS copy: ENA FASTQ (_1/_2 interleaved, or the single file); gzip -t catches truncation
    local ftp f; ftp=$(curl -sS --retry 5 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$run&result=read_run&fields=fastq_ftp&format=tsv" | tail -1 | cut -f2)
    local files=($(tr ';' '\n' <<< "$ftp" | grep -E "/${run}(_[12])?\.fastq\.gz$"))
    (( ${#files[@]} )) || { echo "FAIL $run (not on AWS or ENA)"; return 1; }
    for f in "${files[@]}"; do  # resume dropped connections (-C -); gzip -t is the completeness check
      local try
      for try in 1 2 3 4 5 6; do
        gzip -t $T/${f##*/} 2>/dev/null && break
        curl -sS --retry 5 --retry-all-errors -C - -o $T/${f##*/} "http://$f" || true
      done
      gzip -t $T/${f##*/} 2>/dev/null || { rm -f $T/${f##*/}; echo "FAIL $run (ENA download)"; return 1; }
    done
    t_dl=$SECONDS
    local fq=("${files[@]##*/}"); fq=("${fq[@]/#/$T/}")
    if (( ${#fq[@]} == 2 )); then
      paste <(zcat ${fq[0]} | paste - - - -) <(zcat ${fq[1]} | paste - - - -) | tr '\t' '\n'
    else
      zcat ${fq[0]}
    fi | "$CLADEPRINT" subsample --fraction $frac --seed $SEED --counts $T/$run.count | gzip -1 > $T/$run.fq.gz
    read seen kept _ < $T/$run.count
    rm -f $T/$run.count
    [[ ${KEEP_SRA:-0} == 1 ]] || rm -f "${fq[@]}"
  fi
  [[ $seen -eq $spots ]] || { rm -f $T/$run.fq.gz; echo "FAIL $run (stream had $seen of $spots spots)"; return 1; }
  local t_dec=$SECONDS
  # species screen: ~SCREEN_SPOTS spots (50k: enough for the call) spread over the subsample (not the first ones: some uploads
  # are sorted by genome position); share of reads each species wins (MAPQ >= 20, nuclear contigs),
  # in the order of the panel's species list
  local sf; sf=$(awk -v n=$SCREEN_SPOTS -v k=$kept 'BEGIN{f=n/(k>0?k:1); print (f>1?1:f)}')
  zcat $T/$run.fq.gz | "$CLADEPRINT" subsample --fraction $sf --seed 7 \
    | minimap2 -ax sr -t $KMC_THREADS $panel.mmi - 2>/dev/null | samtools view -F 0x900 - \
    | awk -v run=$run -v nuc=$panel.nuclear.txt -v spf=$panel.species.txt 'BEGIN{OFS="\t"; while((getline l < nuc)>0) ok[l]=1;
        while((getline l < spf)>0) sp[++ns]=l}
      {tot++; if(and($2,4)){un++; next} if($5<20 || !($3 in ok)){lq++; next} split($3,a,"|"); c[a[1]]++}
      END{if(!tot) tot=1; printf "%s\t%d\t%.4f\t%.4f", run, tot, un/tot, lq/tot; for(i=1;i<=ns;i++) printf "\t%.4f", c[sp[i]]/tot; print ""}' \
    > $OUT/species/$run.tsv
  local t_scr=$SECONDS
  "$CLADEPRINT" --threads $KMC_THREADS count --reads $T/$run.fq.gz -o $db --tmp $T/$run.kmc --memory-gb 4 > $OUT/logs/$run.count.txt
  rm -rf $T/$run.fq.gz $T/$run.kmc
  echo -e "$run\t$((t_dl - t0))\t$((t_dec - t_dl))\t$((t_scr - t_dec))\t$((SECONDS - t_scr))" >> $OUT/logs/timing.tsv
  echo "done $run: $kept of $seen spots; $(paste -sd' ' $OUT/logs/$run.count.txt | cut -c1-80)"
}
export -f one; export OUT NSPOTS SEED SCREEN_SPOTS KMC_THREADS PANEL_DIR
[[ ${BASH_SOURCE[0]} == "$0" ]] || return 0

[[ -s $OUT/logs/timing.tsv ]] || echo -e "run\tdownload_s\tdecode_s\tscreen_s\tcount_s" > $OUT/logs/timing.tsv
tail -n +2 ${SHEET:-census_runs.tsv} | awk -F'\t' -v set="${SET:-}" 'set == "" || $3 == set {print $2, $3, $11, $12}' | head -n $LIMIT \
  | xargs -P $JOBS -L1 bash -c 'one "$@" || echo "FAIL $1"' _ | tee -a $OUT/logs/fetch.log
echo "finished: $(ls $OUT/db/*.stats.toml 2>/dev/null | wc -l) databases"
