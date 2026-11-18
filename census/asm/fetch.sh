#!/bin/bash
# Download one NCBI assembly genome FASTA (gzipped) into fa/; skip if present.
set -euo pipefail
acc=$1; out=fa/$acc.fa.gz
[ -s "$out" ] && exit 0
tmp=$(mktemp -d tmp_XXXX)
for try in 1 2 3; do
  if curl -sf --max-time 600 -o $tmp/d.zip "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$acc/download?include_annotation_type=GENOME_FASTA"; then
    unzip -qo $tmp/d.zip -d $tmp && f=$(find $tmp -name '*.fna' | head -1) && [ -n "$f" ] && gzip -c "$f" > $out.part && mv $out.part $out && rm -rf $tmp && echo "OK $acc" && exit 0
  fi
  sleep $((try*10))
done
rm -rf $tmp; echo "FAIL $acc"
