#!/usr/bin/env bash
# Census bundle with Ashton 2019 VNIa sub-clades as a fourth level (literature labels only), run
# before the instance is terminated. balance_children makes VNIa markers reach f_in in every
# sub-clade. usage: subclade_rebuild.sh <pid to wait for>
tail --pid=$1 -f /dev/null
ulimit -n $(ulimit -Hn)
cd ~/cladeprint/census/out/subclade && ln -sfn ../db db
/usr/bin/time -v ~/cladeprint/target/release/cladeprint run config.toml > run.log 2>&1
echo "$(date -u) subclade rebuild exit $? $(grep -E 'Elapsed|Maximum resident' run.log | tr -s ' ' | tr '\n' ' ')" >> ~/r4.log
