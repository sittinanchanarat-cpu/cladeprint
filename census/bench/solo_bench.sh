#!/usr/bin/env bash
# R4 resource benchmark: each tool alone on the instance (m7a.4xlarge, 16 vCPU, 64 GB, gp3 500 MB/s),
# 8 threads, same 204 genomes (BCALM2 unitigs of the presence databases; cladeprint counts them
# itself). Waits until every other job is done. Results: bench/solo/*.time; log ~/r4.log
export PATH=$HOME/.local/bin:$PATH MAMBA_ROOT_PREFIX=$HOME/micromamba
eval "$(micromamba shell hook -s bash)"; micromamba activate bench
ulimit -n $(ulimit -Hn)
log() { echo "$(date -u) solo $*" >> ~/r4.log; }
until grep -q "REQUEUE DONE" ~/r4.log && grep -q "FOCUS TREE DONE" ~/r4.log; do sleep 120; done
B=~/cladeprint/census/out/bench; S=$B/solo; mkdir -p $S; cd $B
T() { /usr/bin/time -v -o $S/$1.time "${@:2}"; local rc=$?; log "$1 exit $rc $(grep -E 'Elapsed|Maximum resident' $S/$1.time | tr -s ' ' | tr '\n' ' ')"; return $rc; }
log "start; load $(cut -d' ' -f1-3 /proc/loadavg); mem avail $(free -g | awk 'NR==2{print $7}') GB"
# 1. cladeprint from FASTA (counting included)
mkdir -p $S/cp && cd $S/cp
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "sample","assembly","species","lineage","sublineage","train"; next} {print $1,"../../genomes/"$1".fa",$3,$4,$5,($6=="train"?"yes":"no")}' $B/samples.tsv > samples.tsv
sed 's/threads = 8/threads = 8/' $B/cp/config.toml > config.toml
T cladeprint ~/cladeprint/target/release/cladeprint run config.toml > run.log 2>&1
cd $B
# 2. Kover: matrix, then SCM (bound) per target
D="sudo docker run --rm -v $B:/data -w /data aldro61/kover:latest"
mkdir -p kover/solo kover/tmp
T kover_create $D kover dataset create from-contigs --genomic-data kover/genomes.tsv --phenotype-description C.neoformans_VNII \
  --phenotype-metadata kover/C.neoformans_VNII.pheno --output kover/solo/base.kover --kmer-size 31 --n-cpu 8 --compression 0 \
  --temp-dir /data/kover/tmp > $S/kover_create.log 2>&1
for t in $(tail -n +2 targets.tsv | cut -f2); do
  sudo cp kover/solo/base.kover kover/solo/$t.kover
  $D python /data/kover_swap.py kover/solo/$t.kover kover/$t.pheno $t > /dev/null 2>&1
  $D kover dataset split --dataset kover/solo/$t.kover --id s --train-ids kover/train.ids --test-ids kover/test.ids --random-seed 1 > /dev/null 2>&1
  T kover_learn_$t $D kover learn scm --dataset kover/solo/$t.kover --split s --model-type conjunction disjunction \
    --p 0.1 0.316 0.562 1.0 1.778 3.162 10.0 999999.0 --hp-choice bound --max-rules 10 --max-equiv-rules 10000 --n-cpu 8 \
    --random-seed 1 --output-dir kover/solo/res_$t > /dev/null 2>&1
  sudo rm -f kover/solo/$t.kover
done
sudo rm -f kover/solo/base.kover
# 3. Neptune per target (defaults, 8 processes); also its R4 output
for t in $(tail -n +2 targets.tsv | cut -f2); do
  [ -s neptune/out_$t/consolidated/consolidated.fasta ] && continue
  T neptune_$t neptune -i $(cat neptune/$t.inc) -e $(cat neptune/$t.exc) -o neptune/out_$t -p 8 > neptune/$t.log 2>&1
  du -sh neptune/out_$t >> ~/r4.log
  find neptune/out_$t -mindepth 1 -maxdepth 1 ! -name consolidated -exec rm -rf {} +
done
log "SOLO DONE"
