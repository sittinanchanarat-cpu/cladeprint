#!/usr/bin/env bash
# Kover matrix build again, alone, for its peak memory (cgroup sampler, ~/docker_peaks.tsv): the first
# solo run's /usr/bin/time measured only the docker client. Runs after the clean cladeprint rerun.
# usage: solo_kover_create_rerun.sh <pid to wait for>
tail --pid=$1 -f /dev/null
B=~/cladeprint/census/out/bench
s=$(date -u +%s)
sudo docker run --rm -v $B:/data -w /data aldro61/kover:latest kover dataset create from-contigs \
  --genomic-data kover/genomes.tsv --phenotype-description C.neoformans_VNII \
  --phenotype-metadata kover/C.neoformans_VNII.pheno --output kover/solo/rerun.kover \
  --kmer-size 31 --n-cpu 8 --compression 0 --temp-dir /data/kover/tmp > $B/solo/kover_create_rerun.log 2>&1
rc=$?
echo "$(date -u) solo kover_create_rerun exit $rc wall $(( $(date -u +%s) - s )) s; peak = last row of ~/docker_peaks.tsv" >> ~/r4.log
sudo rm -f $B/kover/solo/rerun.kover
