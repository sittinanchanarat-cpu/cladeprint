#!/usr/bin/env bash
# R4 tool runs on the benchmark subset (census/bench_r4.py prepare); log ~/r4.log
export PATH=$HOME/.local/bin:$PATH MAMBA_ROOT_PREFIX=$HOME/micromamba
eval "$(micromamba shell hook -s bash)"; micromamba activate bench
ulimit -n $(ulimit -Hn)
B=~/cladeprint/census/out/bench; CP=~/cladeprint/target/release/cladeprint
log() { echo "$(date -u) $*" >> ~/r4.log; }
cd $B
for d in cp cpx; do
  ( cd $d && ln -sfn ../../db db
    /usr/bin/time -v $CP run config.toml > run.log 2>&1 && \
    $CP classify bundle --samples ../test.tsv --kmc-db db --long --tmp tmp_cls -o test_calls.tsv > classify.log 2>&1 )
  log "cladeprint $d exit $?"
done
until grep -q "^EXIT" ~/bench_genomes.log 2>/dev/null; do sleep 60; done
log "genomes: $(tail -1 ~/bench_genomes.log)"
mkdir -p kover/tmp
/usr/bin/time -v -o kover/create.time sudo docker run --rm -v $B:/data -w /data aldro61/kover:latest \
  kover dataset create from-contigs --genomic-data kover/genomes.tsv --output kover/all.kover \
  --kmer-size 31 --n-cpu 8 --compression 0 --temp-dir /data/kover/tmp > kover/create.log 2>&1
log "kover create exit $?"
tail -n +2 targets.tsv | cut -f2 | while read t; do
  /usr/bin/time -v -o neptune/$t.time neptune -i $(cat neptune/$t.inc) -e $(cat neptune/$t.exc) \
    -o neptune/out_$t -p 8 > neptune/$t.log 2>&1
  rc=$?
  log "neptune $t exit $rc $(grep Elapsed neptune/$t.time)"
  du -sh neptune/out_$t >> ~/r4.log
  find neptune/out_$t -mindepth 1 -maxdepth 1 -type d ! -name consolidated -exec rm -rf {} +
done
log "R4 TOOLS DONE"
