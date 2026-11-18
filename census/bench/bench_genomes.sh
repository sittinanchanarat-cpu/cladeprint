#!/usr/bin/env bash
# R4: each benchmark run's presence database -> unitigs (BCALM2), the genome Neptune and Kover see.
# usage: bench_genomes.sh <sample>   (run from census/out; writes bench/genomes/<sample>.fa)
set -euo pipefail
s=$1; d=bench/genomes; w=bench/work/$s
mkdir -p $d $w
[ -s $d/$s.fa ] && exit 0
kmc_tools -t2 transform db/$s dump -s $w/dump.txt
awk '{print ">\n"$1}' $w/dump.txt > $w/kmers.fa
rm $w/dump.txt
( cd $w && bcalm -in kmers.fa -kmer-size 31 -abundance-min 1 -nb-cores 2 -out u > bcalm.log 2>&1 )
awk -v s=$s '/^>/{n++; print ">"s"_u"n; next} {print}' $w/u.unitigs.fa > $d/$s.fa.part
mv $d/$s.fa.part $d/$s.fa
rm -rf $w
