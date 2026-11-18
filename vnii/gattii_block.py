"""VNII markers absent from H99: a ~23 kb C. gattii-type block on the chr7-syntenic locus.

Steps (inputs from kmers/; references from pilot/ref; WM276 annotation fetched once into vnii/):
  1. Unplaced VNII markers: edit-filtered VNII unitigs (kmers/groups/VNII.ed2) with < 80% aligned
     to H99 (same rule as place_markers.py); keep >= 60 bp. Baseline: random H99 nuclear segments
     with the same lengths (2 per unitig).
  2. Map both to the six other species of the complex (minimap2 sr, k 11). Best hit per segment,
     >= 80% aligned; compare identity (VNII-only vs baseline).
  3. Block = VNII segments on WM276 chrF 900-960 kb. Identity of the block to each species separately.
  4. Synteny: the WM276 900-970 kb region aligned to H99 (asm20).
  5. Genes in the block (WM276 GFF).
  6. Presence of the block's k-mers in every isolate (kmers/tmp/present, from score_isolates.py).
Writes vnii/block.fa, vnii/block_kmers.npy, vnii/block_presence.tsv and prints the numbers in README.md.

Usage: micromamba run -n cladeprint python vnii/gattii_block.py
"""
import io
import os
import random
import subprocess
import numpy as np
import pandas as pd

os.chdir(os.path.dirname(os.path.abspath(__file__)))
K = 31
MASK = (1 << 2 * K) - 1
PANEL = "../pilot/ref/panel"
SPECIES = ["gattii", "deuterogattii", "bacillisporus", "tetragattii", "decagattii", "deneoformans"]
CHRF = "gattii|NC_014943.1"
PAF = ["q", "qlen", "qs", "qe", "st", "t", "tlen", "ts", "te", "m", "al", "mq"]


def sh(*cmd):
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout


def paf(text):
    d = pd.read_csv(io.StringIO(text), sep="\t", header=None, usecols=range(12), names=PAF) if text else \
        pd.DataFrame(columns=PAF)
    d = d.sort_values("m", ascending=False).drop_duplicates("q")
    d = d[(d.qe - d.qs) >= 0.8 * d.qlen].copy()
    d["id"] = d.m / d.al
    d["sp"] = d.t.str.split("|").str[0]
    return d


def fasta(path):
    L = open(path).read().split()
    return dict(zip([h[1:] for h in L[::2]], L[1::2]))


def write(path, recs):
    with open(path, "w") as f:
        for h, s in recs:
            f.write(f">{h}\n{s}\n")


def mm(ref, query, *opts):
    return sh("minimap2", "-c", *opts, "--secondary=no", "-t", "3", ref, query)


# 1. unplaced VNII segments and a length-matched H99 baseline
p = pd.read_csv("../kmers/groups/VNII.ed2.paf", sep="\t", header=None, usecols=range(4), names=PAF[:4])
placed = set(p[(p.qe - p.qs) >= 0.8 * p.qlen].q)
uni = fasta("../kmers/groups/VNII.ed2.unitigs.fa")
un = [(h, s) for h, s in uni.items() if h not in placed and len(s) >= 60]
write("unplaced.fa", un)
h99 = {}
for line in open("../pilot/ref/H99.fna"):
    if line[0] == ">":
        name = line[1:].split()[0]
        h99[name] = []
    else:
        h99[name].append(line.strip().upper())
h99 = {k: "".join(v) for k, v in h99.items() if k != "NC_018792.1"}
rng = random.Random(1)
names = list(h99)
w = np.cumsum([len(h99[n]) for n in names]) / sum(len(h99[n]) for n in names)
base = []
for _, s in un:
    for _ in range(2):
        n = names[np.searchsorted(w, rng.random())]
        q = rng.randrange(len(h99[n]) - len(s))
        seg = h99[n][q:q + len(s)]
        if "N" not in seg:
            base.append((f"b{len(base)}", seg))
write("baseline.fa", base)

# 2. VNII-only vs baseline identity to the other species
if not os.path.exists("others.fna"):
    with open("others.fna", "w") as f:
        for sp in SPECIES:
            f.write(open(f"{PANEL}/{sp}.fna").read())
u = paf(mm("others.fna", "unplaced.fa", "-x", "sr", "-k", "11", "-w", "3"))
b = paf(mm("others.fna", "baseline.fa", "-x", "sr", "-k", "11", "-w", "3"))
print(f"VNII-only segments >= 60 bp: {len(un)}; aligned to another species: {len(u)}; baseline {len(b)}/{len(base)}")
for name, d in (("VNII-only", u), ("H99 baseline", b)):
    print(f"  {name}: median identity {d.id.median():.3f}, >= 98%: {(d.id >= 0.98).mean():.0%}; "
          f"C. gattii hits {int((d.sp == 'gattii').sum())} at median {d[d.sp == 'gattii'].id.median():.3f}")
hi = u[u.id >= 0.98]
print(f"  >= 98% hits per target: {hi.t.value_counts().head(3).to_dict()}")

# 3. the block and its identity to each species
blk = u[(u.t == CHRF) & (u.ts > 900_000) & (u.te < 960_000)]
write("block.fa", [(q, uni[q]) for q in blk.q])
print(f"\nBlock: {len(blk)} unitigs, {blk.qlen.sum():,} bp of marker sequence, WM276 chrF "
      f"{blk.ts.min():,}-{blk.te.max():,}")
for sp in SPECIES:
    d = paf(mm(f"{PANEL}/{sp}.fna", "block.fa", "-x", "sr", "-k", "11", "-w", "3"))
    print(f"  {sp:14s} aligned {len(d)}/{len(blk)}, identity {d.m.sum() / max(1, d.al.sum()):.4f}")

# 4. synteny of the WM276 region with H99
reg = sh("samtools", "faidx", f"{PANEL}/gattii.fna", f"{CHRF}:900000-970000").replace(f">{CHRF}:900000-970000", ">reg")
open("wm276_region.fa", "w").write(reg)
print("\nWM276 chrF 900-970 kb vs H99 (asm20):")
for line in mm("../pilot/ref/H99.fna", "wm276_region.fa", "-x", "asm20").splitlines():
    f = line.split("\t")
    print(f"  WM276 {int(f[2]) + 900000:,}-{int(f[3]) + 900000:,} -> H99 {f[5]} {int(f[7]):,}-{int(f[8]):,}"
          f"  id {int(f[9]) / int(f[10]):.3f}")

# 5. genes in the block
if not os.path.exists("wm276.gff"):
    url = ("https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_000185945.1/download"
           "?include_annotation_type=GENOME_GFF")
    subprocess.run(f"curl -sSL --retry 5 '{url}' -o wm276.zip && unzip -p wm276.zip '*.gff' > wm276.gff "
                   "&& rm wm276.zip", shell=True, check=True)
g = pd.read_csv("wm276.gff", sep="\t", comment="#", header=None, usecols=[0, 2, 3, 4, 8])
g = g[(g[0] == CHRF.split("|")[1]) & (g[2] == "CDS") & (g[4] > blk.ts.min()) & (g[3] < blk.te.max())]
g["locus"] = g[8].str.extract(r"locus_tag=([^;]+)")
g["product"] = g[8].str.extract(r"product=([^;]+)")[0].str.replace("%2C", ",")
print("\nGenes overlapping the block:")
print(g.groupby("locus").agg(start=(3, "min"), end=(4, "max"), product=("product", "first")).to_string())


# 6. presence of the block in every isolate
def canon(s):
    out = []
    for i in range(len(s) - K + 1):
        f = 0
        for c in s[i:i + K]:
            f = (f << 2) | "ACGT".index(c)
        r, x = 0, f ^ MASK
        for _ in range(K):
            r, x = (r << 2) | (x & 3), x >> 2
        out.append(min(f, r))
    return out


B = np.unique(np.array([k for q in blk.q for k in canon(uni[q])], np.uint64))
np.save("block_kmers.npy", B)
s = pd.read_csv("../kmers/samples.tsv", sep="\t").set_index("run")
lin = pd.read_csv("../pilot/pilot_lineages.tsv", sep="\t").set_index("run")
pres = pd.DataFrame([dict(run=r, lineage=s.called[r], use=s.use[r], peak_cov=s.peak_cov[r],
                          strain=lin.strain.get(r), country=lin.country.get(r),
                          block=round(float(np.isin(B, np.load(f"../kmers/tmp/present/{r}.npy")).mean()), 4))
                     for r in s.index])
pres.to_csv("block_presence.tsv", sep="\t", index=False)
print(f"\nBlock k-mers: {len(B):,}; share present per isolate:")
print(pres.groupby("lineage").block.describe()[["count", "min", "50%", "max"]].round(3).to_string())
v = pres[(pres.lineage == "VNII") & (pres.use == "train")]
print(f"VNII training isolates: {len(v)} from {v.country.nunique()} countries, block {v.block.min():.3f}-{v.block.max():.3f}")
