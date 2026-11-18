# Shared SRA helpers, sourced by external/run_external.sh and vnii/donor/fetch_vgi.sh (exported for
# their xargs workers). Needs curl, fastq-dump (sra-tools), gzip and the cladeprint binary.

CLADEPRINT=${CLADEPRINT:-$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)/target/release/cladeprint}
DECODE_JOBS=${DECODE_JOBS:-4}  # parallel fastq-dump per run (4: 3.2x, 8: 3.8x faster than one, 15 Oct benchmark)

fetch() {  # url size out: SEGS parallel byte ranges (S3 throttles per connection; 8 ranges ~5x faster)
  local url=$1 size=$2 out=$3 n=${SEGS:-8} i
  local seg=$(( (size + n - 1) / n ))
  for ((i = 0; i < n; i++)); do
    (
      a=$((i * seg)); b=$((a + seg - 1)); (( b >= size )) && b=$((size - 1))
      for try in 1 2 3 4 5 6; do  # resume each segment by appending from where it stopped
        have=$(stat -c %s $out.part$i 2>/dev/null || echo 0)
        (( have >= b - a + 1 )) && break
        curl -sS --retry 5 --retry-all-errors -r $((a + have))-$b "$url" >> $out.part$i || true
      done
    ) &
  done
  wait
  # join only when every part is whole; a short part (network down) is kept for the next fetch
  for ((i = 0; i < n; i++)); do
    a=$((i * seg)); b=$((a + seg - 1)); (( b >= size )) && b=$((size - 1))
    (( $(stat -c %s $out.part$i 2>/dev/null || echo 0) == b - a + 1 )) || return 1
  done
  # join in place: append each part to the first and delete it at once, so the join needs one
  # part of extra space, not a second copy of the run (a 7 GB run overflowed the disk that way)
  mv $out.part0 $out
  for ((i = 1; i < n; i++)); do cat $out.part$i >> $out && rm -f $out.part$i || return 1; done
}

decode() {  # sra spots fraction seed out.fq.gz: subsampled FASTQ; prints "spots_seen spots_kept"
  # The run is decoded as DECODE_JOBS spot ranges in parallel (fastq-dump -N/-X; fasterq-dump is
  # not used: even with -Z it writes the whole FASTQ to its temp directory). `cladeprint subsample`
  # keeps a spot by a hash of its name and the seed, so the sample does not depend on the split;
  # mates share the spot name and stay together. Chunk outputs are separate gzip members.
  local sra=$1 spots=$2 frac=$3 seed=$4 out=$5 p=$DECODE_JOBS i a b
  (( p > spots )) && p=1
  for ((i = 0; i < p; i++)); do
    a=$(( i * spots / p + 1 )); b=$(( (i + 1) * spots / p ))
    ( fastq-dump -N $a -X $b --split-spot --stdout $sra 2>>$out.fqd.log \
        | "$CLADEPRINT" subsample --fraction $frac --seed $seed --counts $out.$i.count | gzip -1 > $out.$i ) &
  done
  wait
  : > $out
  for ((i = 0; i < p; i++)); do cat $out.$i >> $out && rm -f $out.$i; done
  cat $out.*.count | awk '{s += $1; k += $2} END {print s + 0, k + 0}'
  rm -f $out.*.count $out.fqd.log
}
wait_net() {  # url: block until the URL answers (offline, or just resumed from suspend)
  until curl -sSI --max-time 20 "$1" > /dev/null 2>&1; do sleep 60; done
}

export -f fetch decode wait_net; export CLADEPRINT DECODE_JOBS
