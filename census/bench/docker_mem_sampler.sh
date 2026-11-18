#!/usr/bin/env bash
# peak memory of each docker container (cgroup v2 memory.current, polled every 2 s); /usr/bin/time
# around `docker run` only measures the client. Writes $1 (default /root/docker_peaks.tsv): time, id, peak memory.current (incl. page cache), peak anon bytes
declare -A peak peakanon
OUT=${1:-/root/docker_peaks.tsv}
while true; do
  for d in /sys/fs/cgroup/system.slice/docker-*.scope; do
    [ -e "$d/memory.current" ] || continue
    id=$(basename $d .scope); id=${id#docker-}; v=$(cat $d/memory.current 2>/dev/null) || continue
    a=$(awk '$1=="anon"{print $2}' $d/memory.stat 2>/dev/null)
    if [ -z "${peak[$id]}" ] || [ "$v" -gt "${peak[$id]}" ]; then peak[$id]=$v; fi
    if [ -n "$a" ] && { [ -z "${peakanon[$id]}" ] || [ "$a" -gt "${peakanon[$id]}" ]; }; then peakanon[$id]=$a; fi
  done
  for id in "${!peak[@]}"; do
    if [ ! -e /sys/fs/cgroup/system.slice/docker-$id.scope ]; then
      echo -e "$(date -u +%FT%T)\t$id\t${peak[$id]}\t${peakanon[$id]}" >> $OUT; unset peak[$id] peakanon[$id]
    fi
  done
  sleep 2
done
