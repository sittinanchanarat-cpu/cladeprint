#!/usr/bin/env bash
# one R2/R4 job: full run in census/out/r2/<tag>, then classify the held-out runs against it
tag=$1
ulimit -n $(ulimit -Hn)
cd ~/cladeprint/census/out/r2/$tag || exit 1
ln -sfn ../../db db
echo "$(date -u) START $tag" >> ~/r2.log
/usr/bin/time -v ~/cladeprint/target/release/cladeprint run config.toml > run.log 2>&1
rc=$?
if [ $rc -eq 0 ] && [ -f heldout.tsv ]; then
  ~/cladeprint/target/release/cladeprint classify bundle --samples heldout.tsv --kmc-db db --long --tmp tmp_cls -o heldout_calls.tsv > classify.log 2>&1
  rc=$?
fi
echo "$(date -u) END $tag exit $rc" >> ~/r2.log
