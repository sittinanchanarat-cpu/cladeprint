#!/usr/bin/env bash
# Faster variant of align.sh for the CSF runs: a 10M-spot random subsample (scripts/sra.sh decode:
# parallel fastq-dump + `cladeprint subsample`, deterministic per spot), and no coordinate sort -
# paired mates are split into two files; hisat2 output goes to an unsorted BAM (mates stay adjacent) that featureCounts reads directly.
# YPD (SRR1015479) was aligned at full depth by align.sh; analyze.py compares it with a subsample.
# Usage: bash align_fast.sh [RUN...]   (default: both CSF runs, then a YPD subsample to check the depth)
set -euo pipefail
cd "$(dirname "$0")"
source ../../scripts/sra.sh
export PATH=$HOME/micromamba/envs/cladeprint/bin:$PATH  # fastq-dump for decode()
SPOTS_WANTED=10000000
declare -A SPOTS=([SRR1015479]=41283760 [SRR1015480]=26051667 [SRR1015481]=41270151)
mkdir -p sub
runs=("$@"); (( ${#runs[@]} )) || runs=(SRR1015481 SRR1015480 SRR1015479)
for r in "${runs[@]}"; do
  if [[ ! -s sub/$r.counts.tsv ]]; then
    frac=$(awk -v w=$SPOTS_WANTED -v n=${SPOTS[$r]} 'BEGIN {f = w / n; print (f > 1 ? 1 : f)}')
    [[ -s sub/$r.fq.gz ]] || DECODE_JOBS=4 decode sra/$r.sra ${SPOTS[$r]} $frac 20261016 sub/$r.fq.gz > sub/$r.decode.txt
    if [[ $r == SRR1015480 ]]; then  # mates into two files (hisat2 2.2.3 has no --interleaved;
      # named pipes deadlock: hisat2 reads one mate file in large blocks while the splitter waits on the other)
      [[ -s sub/$r.2.fq.gz ]] || zcat sub/$r.fq.gz | paste - - - - - - - - |
        awk -F'\t' -v a="gzip -1 > sub/$r.1.fq.gz" -v b="gzip -1 > sub/$r.2.fq.gz" \
          '{print $1"\n"$2"\n"$3"\n"$4 | a; print $5"\n"$6"\n"$7"\n"$8 | b}'
      mode=(-1 sub/$r.1.fq.gz -2 sub/$r.2.fq.gz); pe=(-p --countReadPairs)
    else mode=(-U sub/$r.fq.gz); pe=(); fi
    micromamba run -n rnaseq hisat2 -p 7 -x vnii_idx "${mode[@]}" --max-intronlen 5000 2> sub/$r.log |
      micromamba run -n rnaseq samtools view -b -o sub/$r.unsorted.bam -
    micromamba run -n rnaseq featureCounts -F SAF -a genes.saf -s 0 -Q 10 -T 4 "${pe[@]}" \
      -o sub/$r.counts.tsv sub/$r.unsorted.bam > sub/$r.fc.log 2>&1
    rm -f sub/$r.fq.gz sub/$r.[12].fq.gz
  fi
  echo "$r: spots seen/kept $(cat sub/$r.decode.txt); $(tail -1 sub/$r.log)"
done
