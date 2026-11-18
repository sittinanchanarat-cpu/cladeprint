"""Fig. 6  VNIa is several clades, and markers for each recover them.

A  neighbour-joining tree of 524 VNI runs (split k-mer distances, census/out/r3/nj.nwk), tips coloured by
   Ashton 2019 VNIa sub-clade, else by this study's nearest-neighbour assignment (census/r3_fig.py logic;
   the two-strip version with candidate runs is the supplementary figure).
B  VNIa marker sensitivity (leave-one-out) of each sub-clade's runs: census bundle (VNIa markers in >= 90% of
   VNIa) vs the sub-clade bundle (also in >= 90% of every sub-clade).
C  leave-one-out sensitivity of each sub-clade's own markers (sub-clade bundle).

Style as make_figures.py (SC Lab). Usage: micromamba run -n cladeprint python figures/fig6_vnia.py
"""
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
SC_ANNOT, SC_AXIS = "#4A4A4A", "#2C2C2C"
mpl.rcParams.update({
    "font.family": ["Arial", "Liberation Sans", "DejaVu Sans"], "font.size": 7, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.75, "axes.edgecolor": SC_AXIS,
    "xtick.direction": "out", "ytick.direction": "out", "xtick.major.size": 3.5, "ytick.major.size": 3.5,
    "xtick.major.width": 0.75, "ytick.major.width": 0.75, "xtick.color": SC_AXIS, "ytick.color": SC_AXIS,
    "lines.linewidth": 1.5, "lines.markersize": 5, "legend.frameon": False, "figure.dpi": 300,
    "savefig.facecolor": "white", "savefig.bbox": "tight", "pdf.fonttype": 42, "ps.fonttype": 42,
})
SUBS = ["VNIa-4", "VNIa-5", "VNIa-93", "VNIa-32", "VNIa-X", "VNIa-Y"]
COLOR = dict(zip(SUBS, ["#9B2335", "#C96A3A", "#DFB166", "#B8C4A2", "#5E7A52", "#2F4A3E"]))
COLOR.update({"VNIa (unassigned)": "#D0D0D0", "VNIb": "#4A4A4A", "VNIc": "#9A9A9A"})


def panel(ax, letter, dx, dy=1.02):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")


# labels as in census/r3_fig.py
L = pd.read_csv("labels/literature_labels.tsv", sep="\t", dtype=str)
L["sub"] = L.label_raw.str.extract(r"(VNIa-[0-9A-Z]+)", expand=False)
lab = L[L["sub"].notna()].assign(run=L.runs.str.split(";")).explode("run").drop_duplicates("run").set_index("run")["sub"]
taxa = pd.read_csv("census/out/r3/taxa.tsv", sep="\t", dtype=str).set_index("sample")
nn = pd.read_csv("census/out/r3/vnia_nn_subclade.tsv", sep="\t", index_col=0)
ambiguous = set(nn.index[(nn.d1 > 300) | (nn.d2 < 3 * nn.d1)])


def group(s):
    sub = taxa.label_prop.get(s)
    if sub in ("VNIb", "VNIc"):
        return sub
    if s in lab.index:
        return lab[s]
    if s in nn.index and s not in ambiguous:
        return nn.loc[s, "sub"]
    return "VNIa (unassigned)"


fig = plt.figure(figsize=(7.0, 4.4), constrained_layout=True)
gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.1])
a = fig.add_subplot(gs[:, 0])
b, c = fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 1])

# A: tree, tips as a coloured strip
tree = Phylo.read("census/out/r3/nj.nwk", "newick")
tips = tree.get_terminals()
y = {t: i for i, t in enumerate(tips)}
depth = tree.depths()


def ypos(cl):
    if cl not in y:
        y[cl] = sum(ypos(k) for k in cl.clades) / len(cl.clades)
    return y[cl]


ypos(tree.root)
xmax = max(depth[t] for t in tips)
for cl in tree.find_clades(order="preorder"):
    for k in cl.clades:
        a.plot([depth[cl], depth[cl]], [y[cl], y[k]], color=SC_AXIS, lw=0.25)
        a.plot([depth[cl], depth[k]], [y[k], y[k]], color=SC_AXIS, lw=0.25)
strip, x0 = xmax * 0.05, xmax * 1.02
for t in tips:
    a.add_patch(mpl.patches.Rectangle((x0, y[t] - 0.5), strip, 1, color=COLOR[group(t.name)], lw=0))
a.set_ylim(-2, len(tips) + 2)
a.set_xlim(0, xmax * 1.1)
a.set_yticks([])
a.spines["left"].set_visible(False)
a.set_xlabel("SNPs per Mb (NJ)")
a.legend(handles=[mpl.patches.Patch(color=v, label=k) for k, v in COLOR.items()], loc="lower left",
         fontsize=6, handlelength=0.9)
panel(a, "A", dx=-0.02)

# B: VNIa markers carried by each sub-clade, before and after requiring them in every sub-clade
rng = np.random.default_rng(1)
med = {}
for k, (f, name, n) in enumerate([("census/out/loo.tsv", "VNIa level only", None),
                                  ("census/out/subclade/loo.tsv", "every sub-clade", None)]):
    t = pd.read_csv(f, sep="\t")
    v = t[t.group == "C.neoformans/VNI/VNIa"].assign(sub=lambda d: d["sample"].map(lab)).dropna(subset=["sub"])
    nm = int(v.markers.iloc[0])
    for i, sub in enumerate(SUBS):
        s = v[v["sub"] == sub].sensitivity.values
        x = i + (k - 0.5) * 0.36
        b.scatter(x + rng.uniform(-0.09, 0.09, len(s)), s, s=3, color=COLOR[sub], alpha=0.35 if k == 0 else 0.85,
                  linewidths=0, zorder=2)
        b.plot([x - 0.13, x + 0.13], [np.median(s)] * 2, color=SC_AXIS, lw=1.1, zorder=3)
        med[(sub, k)] = np.median(s)
    b.scatter([], [], s=10, color="#8A8A8A", alpha=0.35 if k == 0 else 0.85, label=f"{name} ({nm} markers)")
b.set_xticks(range(len(SUBS)), SUBS, rotation=30, ha="right")
b.set(ylim=(0, 1.05), ylabel="VNIa markers present\n(leave-one-out)")
b.legend(loc="lower right", handletextpad=0.1, fontsize=6.5)
panel(b, "B", dx=-0.16)

# C: each sub-clade's own markers
t = pd.read_csv("census/out/subclade/loo.tsv", sep="\t")
g = pd.read_csv("census/out/subclade/groups.tsv", sep="\t").set_index("group")
for i, sub in enumerate(SUBS):
    s = t[t.group == f"C.neoformans/VNI/VNIa/{sub}"].sensitivity.values
    c.scatter(i + rng.uniform(-0.12, 0.12, len(s)), s, s=3, color=COLOR[sub], alpha=0.85, linewidths=0)
    c.plot([i - 0.18, i + 0.18], [np.median(s)] * 2, color=SC_AXIS, lw=1.1)
    c.text(i, 0.56, f"{int(g.loc[f'C.neoformans/VNI/VNIa/{sub}', 'markers']):,}", ha="center", fontsize=6,
           color=SC_ANNOT)
c.text(2.5, 0.6, "markers per sub-clade", ha="center", fontsize=6, color=SC_ANNOT)
c.set_xticks(range(len(SUBS)), SUBS, rotation=30, ha="right")
c.set(ylim=(0.5, 1.02), ylabel="Own sub-clade markers\n(leave-one-out)")
panel(c, "C", dx=-0.16)

for ext, dpi in (("pdf", 600), ("png", 300), ("tiff", 300)):
    fig.savefig(f"figures/fig6_vnia.{ext}", dpi=dpi)
print("written: fig6_vnia;", {f"{s}/{k}": round(v, 3) for (s, k), v in med.items() if s in ("VNIa-32", "VNIa-X", "VNIa-Y")})
