#!/usr/bin/env bash
# R4: Kover per target: dataset (its from-contigs needs the phenotype at creation), split, SCM.
# usage: r4_kover.sh <tag> ...; log ~/r4.log
B=~/cladeprint/census/out/bench
K="sudo docker run --rm -v $B:/data -w /data aldro61/kover:latest kover"
for t in "$@"; do
  mkdir -p $B/kover/tmp
  ds=kover/$t.kover
  [ -f $B/$ds ] || /usr/bin/time -v -o $B/kover/$t.create.time $K dataset create from-contigs \
    --genomic-data kover/genomes.tsv --phenotype-description "$t" --phenotype-metadata kover/$t.pheno \
    --output $ds --kmer-size 31 --n-cpu 4 --compression 0 --temp-dir /data/kover/tmp > $B/kover/$t.create.log 2>&1
  echo "$(date -u) kover create $t exit $? $(grep Elapsed $B/kover/$t.create.time | tr -s ' ')" >> ~/r4.log
  $K dataset split --dataset $ds --id s --train-ids kover/train.ids --test-ids kover/test.ids \
    --random-seed 1 > $B/kover/$t.split.log 2>&1
  /usr/bin/time -v -o $B/kover/$t.learn.time $K learn scm --dataset $ds --split s \
    --model-type conjunction disjunction --p 0.1 0.316 0.562 1.0 1.778 3.162 10.0 999999.0 \
    --hp-choice bound --max-rules 10 --max-equiv-rules 10000 --n-cpu 4 --random-seed 1 \
    --output-dir kover/res_$t > $B/kover/$t.learn.log 2>&1
  echo "$(date -u) kover learn $t exit $? $(grep Elapsed $B/kover/$t.learn.time | tr -s ' ')" >> ~/r4.log
  rm -f $B/$ds
done
