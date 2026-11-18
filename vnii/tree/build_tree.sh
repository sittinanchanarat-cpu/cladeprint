#!/usr/bin/env bash
# Gene tree of the VNII gattii-type block (WM276 chrF 923,314-947,300).
# Sequences: VNII read consensus (5 isolates, ../donor/vnii_*.fa); C. gattii VGI read consensus
# (the 12 isolates vgi_dist.py assigned to gattii, ../donor/tmp/blk_RUN.fa) + the WM276 reference;
# the homologous region of the other panel genomes (assemblies, not reads mapped to WM276, to avoid
# reference bias at 3-4% divergence), with C. deneoformans as outgroup.
# Then MAFFT (auto) and IQ-TREE (ModelFinder, 1000 ultrafast bootstraps, SH-aLRT).
# Usage: micromamba run -n cladeprint bash build_tree.sh
set -euo pipefail
cd "$(dirname "$0")"
P=../../pilot/ref/panel
: > seqs.fa
for f in ../donor/vnii_*.fa; do r=$(basename $f .fa | sed 's/vnii_//'); sed "1s/.*/>VNII_$r/" $f >> seqs.fa; done
tail -n +2 ../donor/vgi_dist.tsv | awk -F'\t' '$4=="gattii"{print $1"\t"$2"\t"$3}' | while IFS=$'\t' read run strain country; do
  name=$(echo "VGI_${country}_${strain}_${run}" | tr ' ()' '___')
  sed "1s/.*/>$name/" ../donor/tmp/blk_$run.fa >> seqs.fa
done
sed '1s/.*/>VGI_Australia_WM276_ref/' ../donor/wm276_block.fa >> seqs.fa
for sp in deuterogattii bacillisporus tetragattii decagattii deneoformans; do
  # best contig = most aligned bases; region = span of its alignments; strand of the longest hit
  read ctg s e strand < <(minimap2 -c -x asm20 $P/$sp.fna ../donor/wm276_block.fa 2>/dev/null | awk '
    {a[$6]+=$10; if($10>b[$6]){b[$6]=$10; st[$6]=$5}; if(!($6 in lo)||$8<lo[$6])lo[$6]=$8; if($9>hi[$6])hi[$6]=$9}
    END{for(c in a) if(a[c]>best){best=a[c]; bc=c}; print bc, lo[bc]+1, hi[bc], st[bc]}')
  seq=$(samtools faidx $P/$sp.fna "$ctg:$s-$e" | tail -n +2 | tr -d '\n')
  [[ $strand == "-" ]] && seq=$(echo "$seq" | rev | tr ACGTacgt TGCAtgca)
  echo -e ">${sp}_ref\n$seq" >> seqs.fa
  echo "$sp: $ctg:$s-$e ($strand), $((e - s + 1)) bp" >&2
done
grep -c '>' seqs.fa >&2
mafft --auto --thread 4 seqs.fa > aln.fa 2> mafft.log
iqtree3 -s aln.fa -m MFP -B 1000 -alrt 1000 -T 4 -o deneoformans_ref --prefix block -redo > iqtree.log 2>&1
grep -E 'Best-fit model|Log-likelihood of the tree' block.iqtree >&2
