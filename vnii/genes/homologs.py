"""Best homologous locus of each block protein in each genome of the species panel (miniprot,
spliced protein-to-genome alignment).

Per genome: % aa identity / % of the protein aligned, at the rank-1 locus; '*' marks frameshifts or
in-frame stops there (a pseudogene relic). An ortholog in the complex is typically >= 80% over
most of the protein; a paralog or domain-only hit is ~25-45%. '-' = no alignment.
Also writes homologs.tsv (with loci).
Usage: python homologs.py   (expects block_proteins.faa here)
"""
import re, subprocess
from collections import defaultdict

PANEL = "../../pilot/ref/panel"
SPECIES = ["neoformans", "deneoformans", "gattii", "deuterogattii", "bacillisporus", "tetragattii", "decagattii"]
TAGS = {"XP_003194826.1": "CGB_F4610W", "XP_003194829.1": "CGB_F4660C", "XP_003195014.1": "CGB_F4690W"}

prot, qlen, acc = {}, defaultdict(int), None
for line in open("block_proteins.faa"):
    if line.startswith(">"):
        acc = line[1:].split()[0]
        prot[acc] = TAGS.get(acc) or next(w for w in line.split() if w.startswith("CGB_"))
    else:
        qlen[acc] += len(line.strip())

cell, rows = defaultdict(dict), []
for sp in SPECIES:
    gff = subprocess.run(["miniprot", "--gff", "-t3", "-N", "3", f"{PANEL}/{sp}.fna", "block_proteins.faa"],
                         capture_output=True, text=True, check=True).stdout
    for l in gff.splitlines():
        f = l.split("\t")
        if len(f) < 9 or f[2] != "mRNA" or "Rank=1;" not in f[8]:
            continue
        q, qs, qe = re.search(r"Target=(\S+) (\d+) (\d+)", f[8]).groups()
        ident = float(re.search(r"Identity=([\d.]+)", f[8]).group(1))
        broken = "Frameshift=" in f[8] or "StopCodon=" in f[8]
        cov = (int(qe) - int(qs) + 1) / qlen[q]
        cell[q][sp] = f"{100 * ident:.0f}/{100 * cov:.0f}{'*' if broken else ''}"
        rows.append((prot[q], sp, f[0].split("|")[1], f[3], f[4], f"{ident:.3f}", f"{cov:.2f}", int(broken)))

print(f"{'gene':12s}" + "".join(f"{s[:6]:>9s}" for s in SPECIES))
for q, tag in prot.items():
    print(f"{tag:12s}" + "".join(f"{cell[q].get(s, '-'):>9s}" for s in SPECIES))
print("cells: % aa identity / % of protein aligned at the best locus; * = frameshift or stop")
with open("homologs.tsv", "w") as o:
    o.write("gene\tspecies\tcontig\tstart\tend\tidentity\tcoverage\tbroken\n")
    for r in rows:
        o.write("\t".join(map(str, r)) + "\n")
