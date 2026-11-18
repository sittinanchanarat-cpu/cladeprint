#!/usr/bin/env bash
# R4: Kover for every other target from one k-mer matrix (base.kover, built for VNII):
# copy, rewrite the phenotype (kover_swap.py), split, SCM. log ~/r4.log
B=~/cladeprint/census/out/bench
D="sudo docker run --rm -v $B:/data -w /data aldro61/kover:latest"
for t in $(tail -n +2 $B/targets.tsv | cut -f2 | grep -vx C.neoformans_VNII); do
  [ -s $B/kover/res_$t/results.json ] && continue
  sudo cp $B/kover/base.kover $B/kover/$t.kover
  $D python /data/kover_swap.py kover/$t.kover kover/$t.pheno $t > $B/kover/$t.swap.log 2>&1 || { echo "$(date -u) kover swap $t FAILED" >> ~/r4.log; continue; }
  $D kover dataset split --dataset kover/$t.kover --id s --train-ids kover/train.ids --test-ids kover/test.ids \
    --random-seed 1 > $B/kover/$t.split.log 2>&1
  /usr/bin/time -v -o $B/kover/$t.learn.time $D kover learn scm --dataset kover/$t.kover --split s \
    --model-type conjunction disjunction --p 0.1 0.316 0.562 1.0 1.778 3.162 10.0 999999.0 \
    --hp-choice bound --max-rules 10 --max-equiv-rules 10000 --n-cpu 4 --random-seed 1 \
    --output-dir kover/res_$t > $B/kover/$t.learn.log 2>&1
  echo "$(date -u) kover learn $t exit $? $(cat $B/kover/$t.swap.log) $(grep Elapsed $B/kover/$t.learn.time | tr -s ' ')" >> ~/r4.log
  sudo rm -f $B/kover/$t.kover
done
echo "$(date -u) KOVER DONE" >> ~/r4.log
