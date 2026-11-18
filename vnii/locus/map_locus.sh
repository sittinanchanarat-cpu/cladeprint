#!/usr/bin/env bash
# Remap pilot reads competitively to H99 + C. gattii WM276 (panel genomes) to read the VNII alleles at the
# HOC1-RGD1 locus. Pilot BAMs hold every read of the 1M-pair subsample, mapped or not.
# Keeps the full remapped BAM per run, plus the read pairs touching either locus window (for local assembly).
set -euo pipefail
cd "$(dirname "$0")"
H99_WIN="neoformans|NC_026751.1:980000-1030000"
WM_WIN="gattii|NC_014943.1:900000-970000"
for run in "$@"; do
  [[ -s $run.bam ]] && continue
  samtools collate -u -O ../../pilot/bam/$run.bam tmp.$run | samtools fastq -n -1 tmp.$run.1.fq -2 tmp.$run.2.fq -0 /dev/null -s /dev/null - 2>/dev/null
  minimap2 -t3 -ax sr ref2.mmi tmp.$run.1.fq tmp.$run.2.fq 2>/dev/null | samtools sort -@2 -m 500M -o $run.bam -
  samtools index $run.bam
  # read names touching either window, any MAPQ; then pull both mates from the fastqs
  samtools view $run.bam "$H99_WIN" "$WM_WIN" | cut -f1 | sort -u > tmp.$run.names
  for m in 1 2; do seqkit grep -f tmp.$run.names tmp.$run.$m.fq > $run.locus.$m.fq 2>/dev/null || \
     awk 'NR==FNR{k[$1];next} FNR%4==1{p=(substr($1,2) in k)} p' tmp.$run.names tmp.$run.$m.fq > $run.locus.$m.fq; done
  rm -f tmp.$run.*
  echo "$run done: $(samtools view -c -F 0x904 $run.bam) primary, $(wc -l < $run.locus.1.fq | awk '{print $1/4}') locus pairs"
done
