# Figure S5 of the paper.
"""Fig. S5  The 23-kb block tree with a VNII tip taken from the long-read assembly (GCA_022832995.1).

Same alignment and model as Fig. 5C with one sequence added, so the two are directly comparable:
MAFFT --auto, IQ-TREE 3 (HKY+F+G4), 1,000 ultrafast bootstraps, rooted on C. deneoformans.
Input: figures/msystems/data/block_asm.contree. Style as fig5_vnii.py.
"""
import os
import matplotlib as mpl
import glob as _glob
import matplotlib.font_manager as _fm
for _f in _glob.glob(os.path.join(os.environ.get("CLADEPRINT_FONTS", ""), "*.ttf")):
    _fm.fontManager.addfont(_f)
import matplotlib.pyplot as plt
import pandas as pd
from Bio import Phylo

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
os.makedirs("figures/msystems/out", exist_ok=True)
SC_ANNOT, SC_AXIS = "#4A4A4A", "#2C2C2C"
mpl.rcParams.update({
    "font.family": ["Arial", "Liberation Sans", "DejaVu Sans"], "font.size": 7, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.75, "axes.edgecolor": SC_AXIS,
    "figure.dpi": 300, "savefig.facecolor": "white", "savefig.bbox": "tight", "pdf.fonttype": 42,
    "mathtext.fontset": "custom", "mathtext.it": "Arial:italic", "mathtext.rm": "Arial",
})
GATTII, LOST = "#C96A3A", "#9B2335"
REVD = "figures/msystems/data/"

tree = Phylo.read(REVD + "block_asm.contree", "newick")
tree.root_with_outgroup("deneoformans_ref")
out = max(tree.root.clades, key=lambda k: k.branch_length or 0)
true_len, out.branch_length = out.branch_length, 0.004
tips = tree.get_terminals()
y = {t: i for i, t in enumerate(tips)}
depth = tree.depths()


def ypos(cl):
    if cl not in y:
        y[cl] = sum(ypos(k) for k in cl.clades) / len(cl.clades)
    return y[cl]


ypos(tree.root)
fig, c = plt.subplots(figsize=(3.6, 3.4))
for cl in tree.find_clades(order="preorder"):
    for k in cl.clades:
        c.plot([depth[cl], depth[cl]], [y[cl], y[k]], color=SC_AXIS, lw=0.6)
        c.plot([depth[cl], depth[k]], [y[k], y[k]], color=SC_AXIS, lw=0.6)
        if not k.is_terminal() and k.confidence is not None and k.confidence >= 95 and len(k.get_terminals()) >= 5:
            c.text(depth[k] - 0.0004, y[k] - 0.35, f"{k.confidence:.0f}", fontsize=5, ha="right", va="bottom",
                   color=SC_ANNOT)

VGI_NAME = {r.run: (r.country, r.strain) for r in pd.read_csv("vnii/donor/vgi_samples.tsv", sep="\t").itertuples()}
SHORT_COUNTRY = {"Papua New Guinea": "PNG", "South Africa": "S. Africa"}


def tip_label(n):
    if n == "VNII_GCA_022832995_asm":
        return "VNII GCA_022832995.1 (assembly)", LOST
    if n.startswith("VNII_"):
        return "VNII " + n[5:], LOST
    if n.startswith("VGI_"):
        run = n.rsplit("_", 1)[1]
        if run in VGI_NAME:
            country, strain = VGI_NAME[run]
            return f"VGI {SHORT_COUNTRY.get(country, country)} {strain}".replace("_", " "), GATTII
        return "VGI " + n[4:].replace("_ref", "").replace("_", " "), GATTII
    return "C. " + n.replace("_ref", ""), SC_AXIS


xmax = max(depth[t] for t in tips)
for t in tips:
    lab, colour = tip_label(t.name)
    c.text(depth[t] + xmax * 0.02, y[t], lab, va="center", fontsize=5.5, color=colour,
           style="italic" if lab.startswith("C. ") else "normal",
           fontweight="bold" if "assembly" in lab else "normal")
c.set_xlim(0, xmax * 1.95)
c.set_ylim(len(tips) + 1.6, -1)
c.set_yticks([])
c.set_xticks([])
for s in ("left", "bottom"):
    c.spines[s].set_visible(False)
ob = ypos(out)
c.text(depth[out] * 0.5, ob, "//", ha="center", va="center", fontsize=7, color=SC_AXIS,
       bbox=dict(fc="white", ec="none", pad=0.2))
c.text(depth[out] * 0.5, ob + 0.9, f"{true_len:.3f}", ha="center", va="top", fontsize=5, color=SC_ANNOT)
c.plot([xmax * 0.95, xmax * 0.95 + 0.002], [len(tips) + 0.3] * 2, color=SC_AXIS, lw=0.8)
c.text(xmax * 0.95 + 0.001, len(tips) + 0.5, "0.002 subst./site", ha="center", va="top", fontsize=5.5)
fig.tight_layout()
for ext, dpi in (("pdf", 600), ("png", 300), ("tiff", 300)):
    fig.savefig(f"figures/msystems/out/FigureS5.{ext}", dpi=dpi)
print("written: figures/msystems/out/FigureS5")
