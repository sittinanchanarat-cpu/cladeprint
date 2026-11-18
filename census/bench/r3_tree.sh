#!/usr/bin/env bash
# R3: VNIa phylogeny. SKA2 split-k-mer alignment (k = 31: each split k-mer is a database 31-mer,
# so the BCALM2 unitigs of the presence databases lose nothing) -> IQ-TREE 3, GTR+F+ASC
# (variable sites only), 1000 ultrafast bootstraps, rooted on VNB. log ~/r4.log
export PATH=$HOME/.local/bin:$PATH MAMBA_ROOT_PREFIX=$HOME/micromamba
eval "$(micromamba shell hook -s bash)"; micromamba activate bench
log() { echo "$(date -u) r3 $*" >> ~/r4.log; }
until grep -q "r3 genomes exit" ~/r4.log; do sleep 60; done
cd ~/cladeprint/census/out/r3
tail -n +2 taxa.tsv | cut -f1 | while read s; do f=../bench/genomes/$s.fa; [ -s $f ] && printf "%s\t%s\n" $s $f; done > ska_files.tsv
log "taxa with genomes: $(wc -l < ska_files.tsv)"
# one build of all 530 genomes held every sample's split k-mers at once (OOM at 57 GB): build in
# batches of 40, then merge (memory ~ split k-mers x samples)
rm -f batch_* vnia.skf
split -l 40 -d ska_files.tsv batch_
for b in batch_??; do
  /usr/bin/time -v -o $b.time ska build -k 31 -f $b -o $b --threads 4 > $b.log 2>&1 || { log "ska build $b FAILED"; exit 1; }
done
log "ska batches built: $(ls batch_??.skf | wc -l), peak $(grep -h 'Maximum resident' batch_??.time | awk '{print $NF}' | sort -n | tail -1) kB"
/usr/bin/time -v -o ska_merge.time ska merge -o vnia batch_??.skf > ska_merge.log 2>&1
log "ska merge exit $? $(grep -E 'Elapsed|Maximum resident' ska_merge.time | tr -s ' ' | tr '\n' ' ')"
/usr/bin/time -v -o ska_align.time ska align --min-freq 0.9 --filter no-ambig-or-const --no-gap-only-sites --threads 8 -o aln.fa vnia.skf > ska_align.log 2>&1
log "ska align exit $? sites $(awk '!/^>/{print length($0); exit}' aln.fa)"
roots=$(awk -F'\t' '$2=="VNB root"{print $1}' taxa.tsv | paste -sd,)
run_iq() { /usr/bin/time -v -o iqtree.time iqtree3 -s $1 -m GTR+F+ASC -B 1000 -T 8 -o "$roots" --prefix vnia -redo > iqtree.out 2>&1; }
run_iq aln.fa; rc=$?
if [ $rc -ne 0 ] && [ -s vnia.varsites.phy ]; then log "iqtree: invariant sites under ASC, rerun on varsites"; run_iq vnia.varsites.phy; rc=$?; fi
log "iqtree exit $rc $(grep -E 'Elapsed|Maximum resident' iqtree.time | tr -s ' ' | tr '\n' ' ')"
log "TREE DONE"
