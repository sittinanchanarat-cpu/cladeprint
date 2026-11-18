#!/usr/bin/env bash
# HC1 (VNII) RNA-seq of Chen et al. 2014 (GSE51573) mapped to the VNII long-read genome
# (GCA_022832995.1), to test whether the C. gattii-type block genes are expressed.
# Needs: sra/<run>.sra (downloaded), vnii_idx (hisat2-build of ../junction/vnii_ont.fna), genes.saf
# (miniprot models: the 11 WM276 block/flank proteins + non-overlapping H99 proteins).
# Usage: bash align.sh   (envs: cladeprint for fastq-dump, rnaseq for hisat2/samtools/featureCounts)
set -euo pipefail
cd "$(dirname "$0")"
declare -A COND=([SRR1015479]=YPD [SRR1015480]=exvivo_CSF [SRR1015481]=invivo_CSF)
for r in SRR1015479 SRR1015481 SRR1015480; do
  [[ -s bam/$r.bam ]] && continue
  until [[ -s sra/$r.sra ]] && ! compgen -G "sra/$r.sra.part*" > /dev/null; do sleep 30; done  # download finished
  mkdir -p bam
  if [[ $r == SRR1015480 ]]; then  # paired 2x100: interleaved stream
    micromamba run -n cladeprint fastq-dump --split-spot --stdout sra/$r.sra 2>/dev/null |
      micromamba run -n rnaseq hisat2 -p 6 -x vnii_idx --interleaved - --max-intronlen 5000 2> bam/$r.log |
      micromamba run -n rnaseq samtools sort -@2 -m 600M -o bam/$r.bam -
  else  # single 36 bp
    micromamba run -n cladeprint fastq-dump --stdout sra/$r.sra 2>/dev/null |
      micromamba run -n rnaseq hisat2 -p 6 -x vnii_idx -U - --max-intronlen 5000 2> bam/$r.log |
      micromamba run -n rnaseq samtools sort -@2 -m 600M -o bam/$r.bam -
  fi
  micromamba run -n rnaseq samtools index bam/$r.bam
  echo "$r (${COND[$r]}) done: $(tail -1 bam/$r.log)"
done
# unstranded, unique alignments only (MAPQ >= 10); pairs counted as fragments
for s in 0 1 2; do
  micromamba run -n rnaseq featureCounts -F SAF -a genes.saf -s $s -Q 10 -T 4 -o counts_single_s$s.tsv \
    bam/SRR1015479.bam bam/SRR1015481.bam > /dev/null 2>&1
done
micromamba run -n rnaseq featureCounts -F SAF -a genes.saf -s 0 -Q 10 -T 4 -p --countReadPairs \
  -o counts_paired_s0.tsv bam/SRR1015480.bam > /dev/null 2>&1
echo counts done
