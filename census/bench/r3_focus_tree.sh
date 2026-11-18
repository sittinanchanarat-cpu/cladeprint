#!/usr/bin/env bash
# R3 focused tree (~93 taxa, census/r3_focus.py): SKA2 (k = 31, two build batches + merge) ->
# IQ-TREE 3 GTR+F+ASC, 1000 ultrafast bootstraps, rooted on VNB. Memory capped at 30 GB
# (ulimit -v) so a misestimate fails this job instead of OOM-killing others. log ~/r4.log
export PATH=$HOME/.local/bin:$PATH MAMBA_ROOT_PREFIX=$HOME/micromamba
eval "$(micromamba shell hook -s bash)"; micromamba activate bench
log() { echo "$(date -u) r3focus $*" >> ~/r4.log; }
ulimit -v 31457280
mkdir -p ~/cladeprint/census/out/r3/focus && cd ~/cladeprint/census/out/r3/focus
tail -n +2 ../focus_taxa.tsv | cut -f1 | while read s; do f=../../bench/genomes/$s.fa; [ -s $f ] && printf "%s\t%s\n" $s $f; done > files.tsv
log "taxa with genomes: $(wc -l < files.tsv)"
split -n l/2 -d files.tsv batch_
for b in batch_0?; do
  /usr/bin/time -v -o $b.time ska build -k 31 -f $b -o $b --threads 4 > $b.log 2>&1 || { log "build $b FAILED"; exit 1; }
done
/usr/bin/time -v -o merge.time ska merge -o focus batch_0?.skf > merge.log 2>&1 || { log "merge FAILED $(grep 'Maximum resident' merge.time)"; exit 1; }
log "merge ok $(grep -E 'Elapsed|Maximum resident' merge.time | tr -s ' ' | tr '\n' ' ')"
/usr/bin/time -v -o align.time ska align --min-freq 0.9 --filter no-ambig-or-const --no-gap-only-sites --threads 6 -o aln.fa focus.skf > align.log 2>&1 || { log "align FAILED"; exit 1; }
log "align ok: $(grep -c '>' aln.fa) taxa, $(awk '!/^>/{n+=length($0)} /^>/{if(n){print n; exit}}' aln.fa) sites $(grep 'Maximum resident' align.time | tr -s ' ')"
roots=$(awk -F'\t' '$2=="VNB root"{print $1}' ../focus_taxa.tsv | paste -sd,)
iq() { /usr/bin/time -v -o iqtree.time iqtree3 -s $1 -m GTR+F+ASC -B 1000 -T 6 -o "$roots" --prefix focus -redo > iqtree.out 2>&1; }
iq aln.fa; rc=$?
if [ $rc -ne 0 ] && [ -s focus.varsites.phy ]; then log "invariant sites under ASC: rerun on varsites"; iq focus.varsites.phy; rc=$?; fi
log "iqtree exit $rc $(grep -E 'Elapsed|Maximum resident' iqtree.time | tr -s ' ' | tr '\n' ' ')"
log "FOCUS TREE DONE"
