#!/usr/bin/env bash
# R4: Kover with 5-fold cross-validated hyperparameters (--hp-choice cv), same matrix and
# train/test split as r4_kover2.sh (bound selection). log ~/r4.log
B=~/cladeprint/census/out/bench
D="sudo docker run --rm -v $B:/data -w /data aldro61/kover:latest"
for t in $(tail -n +2 $B/targets.tsv | cut -f2); do
  [ -s $B/kover/res_cv_$t/results.json ] && continue
  sudo cp $B/kover/base.kover $B/kover/$t.kover
  $D python /data/kover_swap.py kover/$t.kover kover/$t.pheno $t > $B/kover/$t.swap.log 2>&1 || { echo "$(date -u) kover swap $t FAILED" >> ~/r4.log; continue; }
  $D kover dataset split --dataset kover/$t.kover --id s --train-ids kover/train.ids --test-ids kover/test.ids \
    --folds 5 --random-seed 1 > $B/kover/$t.split.log 2>&1
  /usr/bin/time -v -o $B/kover/$t.learn.time $D kover learn scm --dataset kover/$t.kover --split s \
    --model-type conjunction disjunction --p 0.1 0.316 0.562 1.0 1.778 3.162 10.0 999999.0 \
    --hp-choice cv --max-rules 10 --max-equiv-rules 10000 --n-cpu 4 --random-seed 1 \
    --output-dir kover/res_cv_$t > $B/kover/$t.learn.log 2>&1
  echo "$(date -u) kover cv $t exit $? $(cat $B/kover/$t.swap.log) $(grep Elapsed $B/kover/$t.learn.time | tr -s ' ')" >> ~/r4.log
  sudo rm -f $B/kover/$t.kover
done
echo "$(date -u) KOVER CV DONE" >> ~/r4.log
