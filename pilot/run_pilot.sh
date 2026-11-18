#!/usr/bin/env bash
# Pilot lineage genotyping: fetch reads -> random 1M-pair subsample -> map to H99 -> joint haploid calls.
#
# Reads come from the SRA copy on AWS open data (free, no account; ~4x faster than ENA from here).
# Whole runs are downloaded and size-checked, then randomly subsampled: taking the first N reads of
# a stream is unsafe because streams get cut off and many uploads are sorted by genome position.
#
# Usage: micromamba run -n cladeprint bash run_pilot.sh            (runs in pilot_samples.tsv)
#        SAMPLES="pilot_samples.tsv ../labels/validation_samples.tsv" micromamba run ...  (several sheets)
#        REDO="SRR1 SRR2" micromamba run -n cladeprint bash run_pilot.sh   (re-fetch these even if a BAM exists)
# Needs ref/ from fetch_refs.sh. Samples with a BAM are skipped; a failed sample keeps its old BAM.
set -euo pipefail
cd "$(dirname "$0")"
REF=ref/H99.fna
NPAIRS=${NPAIRS:-1000000}
JOBS=${JOBS:-2}
REDO=${REDO:-}
SAMPLES=${SAMPLES:-pilot_samples.tsv}
mkdir -p fastq bam logs

fetch_and_map() {
  run=$1
  [[ -s bam/$run.bam && " $REDO " != *" $run "* ]] && return 0
  url="https://sra-pub-run-odp.s3.amazonaws.com/sra/$run/$run"
  size=$(curl -sSI --retry 5 "$url" | tr -d '\r' | awk 'tolower($1)=="content-length:"{print $2}')
  [[ $size =~ ^[0-9]+$ ]] || { echo "FAIL $run (not on AWS)"; return 1; }
  for try in 1 2 3; do
    curl -sS --retry 5 --retry-all-errors -C - -o fastq/$run.sra "$url" || true
    [[ $(stat -c %s fastq/$run.sra 2>/dev/null || echo 0) -eq $size ]] && break
  done
  [[ $(stat -c %s fastq/$run.sra 2>/dev/null || echo 0) -eq $size ]] || { rm -f fastq/$run.sra; echo "FAIL $run (download incomplete)"; return 1; }
  rm -rf fastq/$run.tmp
  fasterq-dump --split-files --threads 2 -t fastq/$run.tmp -O fastq fastq/$run.sra >logs/$run.fqd.log 2>&1 \
    || { rm -rf fastq/$run.sra fastq/$run.tmp fastq/${run}_*.fastq; echo "FAIL $run (fasterq-dump)"; return 1; }
  rm -rf fastq/$run.sra fastq/$run.tmp
  ls -la fastq/${run}* >> logs/$run.fqd.log
  l1=$(wc -l < fastq/${run}_1.fastq); l2=$(wc -l < fastq/${run}_2.fastq)
  [[ $l1 -eq $l2 && $l1 -gt 0 ]] || { rm -f fastq/${run}_*.fastq; echo "FAIL $run (R1 $l1 / R2 $l2 lines)"; return 1; }
  spots=$((l1 / 4))
  frac=$(awk -v n=$NPAIRS -v s=$spots 'BEGIN{f=1.05*n/s; print (f>1?1:f)}')
  if [[ $frac == 1 ]]; then  # fewer spots than NPAIRS: keep all (seqtk reads "1" as a count of 1 read)
    gzip -1 -c fastq/${run}_1.fastq > fastq/${run}_1.fq.gz; gzip -1 -c fastq/${run}_2.fastq > fastq/${run}_2.fq.gz
  else
    seqtk sample -s 20261012 fastq/${run}_1.fastq $frac | gzip -1 > fastq/${run}_1.fq.gz
    seqtk sample -s 20261012 fastq/${run}_2.fastq $frac | gzip -1 > fastq/${run}_2.fq.gz
  fi
  rm -f fastq/${run}_[12].fastq
  minimap2 -ax sr -t 2 -R "@RG\tID:$run\tSM:$run" $REF.mmi fastq/${run}_1.fq.gz fastq/${run}_2.fq.gz 2>logs/$run.mm2.log \
    | samtools sort -@1 -m 400M -o bam/$run.tmp.bam - && samtools index bam/$run.tmp.bam \
    && mv bam/$run.tmp.bam bam/$run.bam && mv bam/$run.tmp.bam.bai bam/$run.bam.bai
  rm -f fastq/${run}_[12].fq.gz
  [[ -s bam/$run.bam ]] || { echo "FAIL $run (no bam)"; return 1; }
  echo "done $run $(samtools view -c -F 0x904 bam/$run.bam) mapped, frac $frac of $spots"
}
export -f fetch_and_map; export REF NPAIRS REDO

for f in $SAMPLES; do tail -n +2 "$f" | cut -f2; done | xargs -P $JOBS -I{} bash -c 'fetch_and_map {} || echo "FAIL {}"' \
  | tee -a logs/fetch_map.log

ls bam/*.bam | grep -v tmp > bam.list
mkdir -p vcf
cut -f1 $REF.fai | xargs -P 4 -I{} bash -c \
  "bcftools mpileup -f $REF -b bam.list -r {} -a AD,DP -q 20 -Q 20 -d 500 -Ou 2>logs/mpileup.{}.log \
   | bcftools call -m --ploidy 1 -v -Oz -o vcf/{}.vcf.gz 2>logs/call.{}.log"
bcftools concat $(cut -f1 $REF.fai | sed 's|^|vcf/|; s|$|.vcf.gz|') -Oz -o calls.raw.vcf.gz.new
mv calls.raw.vcf.gz.new calls.raw.vcf.gz
bcftools index -f calls.raw.vcf.gz
echo "calls done: $(bcftools view -H calls.raw.vcf.gz | wc -l) variant sites"
