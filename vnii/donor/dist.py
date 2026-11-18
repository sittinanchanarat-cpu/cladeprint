"""Divergence of block sequences (VNII consensus per isolate, WM276) from each C. gattii assembly.
minimap2 asm20, all alignments of each query; differences per kb of aligned sequence, where a
difference is a substitution or one indel event, and N positions in the query are skipped.
Usage: python dist.py queries.fa genome.fna[.gz] ...  > dist.tsv"""
import re, subprocess, sys
from collections import defaultdict
q = sys.argv[1]
print("assembly\tquery\taligned_bp\tdiffs_per_kb\tcontigs")
for g in sys.argv[2:]:
    out = subprocess.run(["minimap2", "-c", "--cs", "-x", "asm20", "-t3", g, q], capture_output=True, text=True).stdout
    al, df, ctg = defaultdict(int), defaultdict(int), defaultdict(set)
    for l in out.splitlines():
        f = l.split("\t")
        if "tp:A:P" not in f:
            continue
        cs = next(t for t in f if t.startswith("cs:Z:"))[5:]
        for op, val in re.findall(r"([:*+\-])([0-9a-z]+)", cs):
            if op == ":": al[f[0]] += int(val)
            elif op == "*":
                if val[1] == "n": continue  # N in the query (low coverage in the consensus)
                al[f[0]] += 1; df[f[0]] += 1
            else: df[f[0]] += 1
        ctg[f[0]].add(f[5])
    name = re.sub(r"GCA_\d+\.\d_|\.fna.*|_V1|Cryp_gatt_|Cryp_deca_", "", g.split("/")[-1])
    for k in sorted(al):
        print(f"{name}\t{k}\t{al[k]}\t{1000 * df[k] / al[k]:.1f}\t{len(ctg[k])}")
