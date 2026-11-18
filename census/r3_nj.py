"""R3: neighbour-joining tree of VNI from the census split-k-mer distances (replaces the SKA2 +
IQ-TREE alignment tree, whose merge step needed > 60 GB for 530 taxa).

Distances: census/out/clusters_neo/pairs.tsv (cladeprint clusters; SNPs per Mb of shared covered
genome, exact up to 3,000 per Mb). Taxa: census/out/r3/taxa.tsv minus the VNB rooting runs (beyond
the exact range from VNI) and any run missing more than MAX_MISSING of its distances (partial or
error-rich runs). Remaining gaps are filled with the larger of 3,000 and the run's largest measured
distance (only VNI runs remain, whose pairwise distances are all far below that, so this only
affects excluded-quality runs). Tree: RapidNJ, midpoint-rooted (no outgroup in range).

Checks: each Ashton 2019 sub-clade (literature label) monophyletic or not, and which taxa break it;
the sub-clade clade each unlabelled or candidate run falls in.

Usage: python census/r3_nj.py census/out   -> census/out/r3/nj.nwk, nj_checks.tsv
"""

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo

MAX_MISSING = 0.05


def main(out):
    out = Path(out)
    r3 = out / "r3"
    taxa = pd.read_csv(r3 / "taxa.tsv", sep="\t", dtype=str)
    taxa = taxa[taxa.role != "VNB root"]
    p = pd.read_csv(out / "clusters_neo/pairs.tsv", sep="\t", usecols=[0, 1, 5]).set_axis(["a", "b", "d"], axis=1)
    names = list(taxa["sample"])
    ix = {s: i for i, s in enumerate(names)}
    p = p[p.a.isin(ix) & p.b.isin(ix)]
    n = len(names)
    D = np.full((n, n), np.nan)
    D[p.a.map(ix), p.b.map(ix)] = p.d.values
    D[p.b.map(ix), p.a.map(ix)] = p.d.values
    np.fill_diagonal(D, 0)
    miss = np.isnan(D).mean(axis=1)
    keep = miss <= MAX_MISSING
    dropped = [(names[i], round(miss[i], 3)) for i in np.where(~keep)[0]]
    names = [s for s, k in zip(names, keep) if k]
    D = D[np.ix_(keep, keep)]
    fill = np.maximum(3000, np.nanmax(D, axis=1))
    for i in range(len(names)):
        m = np.isnan(D[i])
        D[i, m] = fill[i]
    D = np.maximum(D, D.T)
    print(f"taxa {len(names)} (dropped {len(dropped)} with > {MAX_MISSING:.0%} missing distances: {dropped})")
    print(f"filled gaps: {int(np.isnan(np.where(np.isnan(D), np.nan, 0)).sum())}")
    with open(r3 / "nj.phy", "w") as f:
        f.write(f"{len(names)}\n")
        for s, row in zip(names, D):
            f.write(s + " " + " ".join(f"{x:.3f}" for x in row) + "\n")
    subprocess.run(["rapidnj", str(r3 / "nj.phy"), "-i", "pd", "-o", "t", "-x", str(r3 / "nj_raw.nwk")],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    tree = Phylo.read(r3 / "nj_raw.nwk", "newick")
    for c in tree.get_terminals():
        c.name = c.name.strip("'")
    tree.root_at_midpoint()
    Phylo.write(tree, r3 / "nj.nwk", "newick")

    # checks against Ashton 2019 sub-clades
    L = pd.read_csv("labels/literature_labels.tsv", sep="\t", dtype=str)
    L["sub"] = L.label_raw.str.extract(r"(VNIa-[0-9A-Z]+)", expand=False)
    lab = L[L["sub"].notna()].assign(run=L.runs.str.split(";")).explode("run").set_index("run")["sub"]
    lab = lab[~lab.index.duplicated()]
    tips = {c.name: c for c in tree.get_terminals()}
    rows = []
    for sub in sorted(lab[lab.index.isin(tips)].unique()):
        members = [tips[s] for s in lab.index[lab == sub] if s in tips]
        mrca = tree.common_ancestor(members) if len(members) > 1 else members[0]
        inside = [t.name for t in mrca.get_terminals()]
        intruders = [s for s in inside if lab.get(s) not in (sub, None) or s not in lab.index]
        other_lab = [s for s in inside if s in lab.index and lab[s] != sub]
        rows.append(dict(subclade=sub, n_labelled=len(members), clade_size=len(inside),
                         monophyletic=not other_lab, other_labelled_inside=",".join(other_lab),
                         n_unlabelled_inside=len([s for s in inside if s not in lab.index])))
    chk = pd.DataFrame(rows)
    chk.to_csv(r3 / "nj_checks.tsv", sep="\t", index=False)
    print(chk.to_string(index=False))
    # where the candidates and outside-scheme runs fall: nearest labelled sub-clade by tree distance
    cand = taxa[taxa.role.str.contains("candidate") | taxa.role.str.contains("small")]["sample"]
    nn = pd.read_csv(r3 / "vnia_nn_subclade.tsv", sep="\t", index_col=0)
    odd = list(dict.fromkeys(list(cand) + [s for s in nn.index[(nn.d1 > 300) | (nn.d2 < 3 * nn.d1)]] + ["ERR842663"]))
    labelled_tips = [t for t in lab.index if t in tips]
    print("\nodd runs: label, nearest labelled tip on the tree (sub-clade, tree distance)")
    for s in odd:
        if s not in tips:
            print(f"  {s}: not in tree")
            continue
        best = min(labelled_tips, key=lambda t: tree.distance(tips[s], tips[t]) if t != s else 1e18)
        print(f"  {s}: {lab.get(s, '-')}, nearest {best} ({lab[best]}, {tree.distance(tips[s], tips[best]):.0f})")


if __name__ == "__main__":
    main(sys.argv[1])
