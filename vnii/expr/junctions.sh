#!/usr/bin/env bash
# Splice support for the block gene models: RNA junctions (>= 3 spliced reads, MAPQ >= 10) in the
# full-depth YPD alignment (bam/SRR1015479.bam) vs the introns of the miniprot models (genes.saf;
# CDS-based, so a junction within 3 bp counts). Writes ypd_block_junctions.txt, block_introns.txt.
set -euo pipefail
cd "$(dirname "$0")"
micromamba run -n rnaseq samtools view -q 10 bam/SRR1015479.bam CP091251.1:404500-428000 |
  awk '$6 ~ /N/ {pos = $4; cig = $6
    while (match(cig, /^[0-9]+[MIDNSH]/)) { n = substr(cig, 1, RLENGTH - 1) + 0; op = substr(cig, RLENGTH, 1)
      if (op == "N") print pos - 1 "\t" pos + n
      if (op == "M" || op == "D" || op == "N") pos += n
      cig = substr(cig, RLENGTH + 1) } }' | sort | uniq -c | awk '$1 >= 3' > ypd_block_junctions.txt
awk -F'\t' '$1 ~ /^BLOCK/ && $1 !~ /4833|5015/' genes.saf | sort -k1,1 -k3,3n |
  awk '{ if ($1 == g) print g "\t" pe "\t" $3 - 1; g = $1; pe = $4 }' > block_introns.txt
python3 - <<'PY'
J = {}
for l in open("ypd_block_junctions.txt"):
    n, a, b = l.split(); J[(int(a), int(b))] = int(n)
tot, hit = {}, {}
for l in open("block_introns.txt"):
    g, a, b = l.split(); a, b = int(a), int(b)
    best = max([n for (x, y), n in J.items() if abs(x - a) <= 3 and abs(y - 1 - b) <= 3] or [0])
    tot[g] = tot.get(g, 0) + 1; hit[g] = hit.get(g, 0) + (best > 0)
print("introns with spliced reads per gene:", {g.split(":")[1]: f"{hit[g]}/{tot[g]}" for g in tot})
print(f"total: {sum(hit.values())}/{sum(tot.values())} model introns supported; {len(J)} junctions "
      f"(>= 3 reads) in the block, {sum(J.values())} spliced reads")
PY
