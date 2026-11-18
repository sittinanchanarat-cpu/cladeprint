"""Is the VNII block one unit? Sliding-window trees and distances along the block alignment.

Windows of 2,000 alignment columns, step 1,000 (aln.fa from build_tree.sh). Per window:
  - mean pairwise p-distance (gaps and N skipped) between groups: VNII vs VGI, VGI within, and
    VNII vs each other species; a mosaic would show windows where VNII is closer to another
    species, or as close to WM276-type VGI as VGI lineages are to each other;
  - IQ-TREE (HKY+F+G4, 1,000 ultrafast bootstraps, rooted on C. deneoformans): is VNII
    monophyletic, and is VNII + VGI a clade with VGI monophyletic (VNII sister to all VGI)?
    Support = ultrafast bootstrap of that clade.
Writes windows.tsv. Usage: micromamba run -n cladeprint python windows.py
"""
import os, re, subprocess
import numpy as np
import pandas as pd

os.chdir(os.path.dirname(os.path.abspath(__file__)))
WIN, STEP = 2000, 1000

names, seqs = [], []
for line in open("aln.fa"):
    if line.startswith(">"):
        names.append(line[1:].strip())
        seqs.append([])
    else:
        seqs[-1].append(line.strip().upper())
aln = np.array([list("".join(s)) for s in seqs])
groups = {n: ("VNII" if n.startswith("VNII") else "VGI" if n.startswith("VGI") else n.replace("_ref", "")) for n in names}
others = ["deuterogattii", "bacillisporus", "tetragattii", "decagattii", "deneoformans"]


def pdist(a, b):
    ok = np.isin(a, list("ACGT")) & np.isin(b, list("ACGT"))
    return (a[ok] != b[ok]).mean() if ok.sum() > 200 else np.nan


def mean_between(cols, g1, g2):
    i1 = [i for i, n in enumerate(names) if groups[n] == g1]
    i2 = [i for i, n in enumerate(names) if groups[n] == g2]
    d = [pdist(aln[a, cols], aln[b, cols]) for a in i1 for b in i2 if a != b]
    return np.nanmean(d) if d else np.nan


def parse(s):
    """Newick -> nested dicts with name (IQ-TREE puts 'SH/UFB' support in internal names)."""
    i = 0

    def node():
        nonlocal i
        ch = []
        if s[i] == "(":
            i += 1
            while True:
                ch.append(node())
                if s[i] == ",":
                    i += 1
                    continue
                i += 1
                break
        m = re.match(r"[^:,();]*", s[i:])
        name = m.group(0)
        i += len(name)
        if s[i] == ":":
            i += len(re.match(r":[0-9.eE-]+", s[i:]).group(0))
        return dict(name=name, ch=ch)

    return node()


def clades(n):
    """(leaf set, support label) for every internal node."""
    if not n["ch"]:
        return {n["name"]}, []
    leaves, out = set(), []
    for c in n["ch"]:
        l, o = clades(c)
        leaves |= l
        out += o
    out.append((frozenset(leaves), n["name"]))
    return leaves, out


vnii = frozenset(n for n in names if groups[n] == "VNII")
vgi = frozenset(n for n in names if groups[n] == "VGI")
rows = []
os.makedirs("win", exist_ok=True)
for start in range(0, aln.shape[1] - WIN + 1, STEP):
    cols = slice(start, start + WIN)
    fa = f"win/w{start}.fa"
    with open(fa, "w") as f:
        for n, row in zip(names, aln[:, cols]):
            f.write(f">{n}\n{''.join(row)}\n")
    subprocess.run(["iqtree3", "-s", fa, "-m", "HKY+F+G4", "-B", "1000", "-T", "2", "-o", "deneoformans_ref",
                    "--prefix", fa[:-3], "-redo", "--quiet"], check=True, capture_output=True)
    _, cl = clades(parse(open(fa[:-3] + ".treefile").read().strip()))
    sup = {c: lab for c, lab in cl}
    ufb = lambda c: (sup[c].split("/")[-1] if c in sup else "-")
    r = dict(start=start, end=start + WIN,
             vnii_mono=vnii in sup, vnii_ufb=ufb(vnii),
             vnii_sister_vgi=(vnii | vgi) in sup and vgi in sup, vnii_vgi_ufb=ufb(vnii | vgi), vgi_ufb=ufb(vgi),
             d_vnii_vgi=mean_between(cols, "VNII", "VGI"), d_vgi_vgi=mean_between(cols, "VGI", "VGI"))
    for o in others:
        r[f"d_vnii_{o}"] = mean_between(cols, "VNII", o)
    r["closest_other_to_vnii"] = min(others, key=lambda o: r[f"d_vnii_{o}"])
    rows.append(r)
    print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items() if not k.startswith("d_vnii_") or k == "d_vnii_vgi"}, flush=True)
df = pd.DataFrame(rows)
df.to_csv("windows.tsv", sep="\t", index=False, float_format="%.4f")
cols = ["start", "end", "vnii_mono", "vnii_ufb", "vnii_sister_vgi", "vnii_vgi_ufb", "vgi_ufb", "d_vnii_vgi", "d_vgi_vgi",
        "d_vnii_tetragattii", "d_vnii_deneoformans"]
print(df[cols].to_string(index=False, float_format="%.4f"))
