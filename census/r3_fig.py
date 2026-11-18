"""R3 supplementary figure: NJ tree of 524 VNI runs (census/r3_nj.py) with two tip strips,
the Ashton 2019 VNIa sub-clade (literature) and this study's nearest-neighbour assignment
(census/out/r3/vnia_nn_subclade.tsv), plus the census mislabel/intermediate candidates.
SC figure style (Desert & Moss palette, Arial, journal double column).

Usage: python census/r3_fig.py census/out   -> figures/r3_vnia_nj.pdf / .tiff
"""

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from Bio import Phylo

SC = ["#9B2335", "#C96A3A", "#DFB166", "#B8C4A2", "#5E7A52", "#2F4A3E"]
SC_AXIS = "#2C2C2C"
mpl.rcParams.update({
    "font.family": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 7, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False, "axes.spines.left": False,
    "axes.linewidth": 0.75, "axes.edgecolor": SC_AXIS, "xtick.direction": "out",
    "xtick.major.size": 3.5, "xtick.major.width": 0.75, "xtick.color": SC_AXIS,
    "legend.frameon": False, "figure.facecolor": "white", "savefig.facecolor": "white",
    "pdf.fonttype": 42, "ps.fonttype": 42,
})
GROUP_COLOR = {
    "VNIa-4": SC[0], "VNIa-5": SC[1], "VNIa-93": SC[2], "VNIa-32": SC[3], "VNIa-X": SC[4],
    "VNIa-Y": SC[5], "VNIa (unassigned)": "#D0D0D0", "VNIb": "#4A4A4A", "VNIc": "#9A9A9A",
}
CANDIDATES = {"SRR822918", "SRR836604", "SRR798180", "SRR1944858", "SRR642856", "ERR2624237",
              "SRR1063253", "ERR842663"}


def main(out):
    out = Path(out)
    tree = Phylo.read(out / "r3/nj.nwk", "newick")
    taxa = pd.read_csv(out / "r3/taxa.tsv", sep="\t", dtype=str).set_index("sample")
    L = pd.read_csv("labels/literature_labels.tsv", sep="\t", dtype=str)
    L["sub"] = L.label_raw.str.extract(r"(VNIa-[0-9A-Z]+)", expand=False)
    lab = L[L["sub"].notna()].assign(run=L.runs.str.split(";")).explode("run").set_index("run")["sub"]
    lab = lab[~lab.index.duplicated()]
    nn = pd.read_csv(out / "r3/vnia_nn_subclade.tsv", sep="\t", index_col=0)
    ambiguous = set(nn.index[(nn.d1 > 300) | (nn.d2 < 3 * nn.d1)])

    def base(s):
        sub = taxa.label_prop.get(s)
        return sub if sub in ("VNIb", "VNIc") else None

    def lit_group(s):
        return lab.get(s) or base(s)

    def our_group(s):
        if base(s):
            return base(s)
        if s in lab.index:
            return lab[s]
        if s in nn.index and s not in ambiguous:
            return nn.loc[s, "sub"]
        return "VNIa (unassigned)"

    tips = tree.get_terminals()
    y = {t: i for i, t in enumerate(tips)}
    depth = tree.depths()

    def ypos(c):
        if c in y:
            return y[c]
        y[c] = sum(ypos(k) for k in c.clades) / len(c.clades)
        return y[c]

    ypos(tree.root)
    xmax = max(depth[t] for t in tips)
    fig, ax = plt.subplots(figsize=(7.0, 9.0))
    for c in tree.find_clades(order="preorder"):
        for k in c.clades:
            ax.plot([depth[c], depth[c]], [y[c], y[k]], color=SC_AXIS, lw=0.3)
            ax.plot([depth[c], depth[k]], [y[k], y[k]], color=SC_AXIS, lw=0.3)
    strip = xmax * 0.03
    x0 = xmax * 1.02
    for t in tips:
        g1, g2 = lit_group(t.name), our_group(t.name)
        if g1:
            ax.add_patch(mpl.patches.Rectangle((x0, y[t] - 0.5), strip, 1, color=GROUP_COLOR[g1], lw=0))
        ax.add_patch(mpl.patches.Rectangle((x0 + 1.3 * strip, y[t] - 0.5), strip, 1, color=GROUP_COLOR[g2], lw=0))
        if t.name in CANDIDATES:
            ax.plot(x0 + 2.6 * strip, y[t], marker="<", color=SC_AXIS, ms=3.5)
            ax.text(x0 + 3.0 * strip, y[t], t.name, va="center", fontsize=5.5)
    ax.text(x0 + 0.5 * strip, len(tips) + 4, "Ashton\n2019", ha="center", va="bottom", fontsize=6)
    ax.text(x0 + 1.8 * strip, len(tips) + 4, "this\nstudy", ha="center", va="bottom", fontsize=6)
    ax.set_ylim(-2, len(tips) + 14)
    ax.set_xlim(0, xmax * 1.25)
    ax.set_yticks([])
    ax.set_xlabel("NJ distance (SNPs per Mb of shared genome)")
    handles = [mpl.patches.Patch(color=c, label=g) for g, c in GROUP_COLOR.items()]
    handles.append(mpl.lines.Line2D([], [], marker="<", ls="", color=SC_AXIS, ms=4, label="census candidate"))
    ax.legend(handles=handles, loc="lower left", ncol=1)
    Path("figures").mkdir(exist_ok=True)
    fig.savefig("figures/r3_vnia_nj.pdf", dpi=600, bbox_inches="tight")
    fig.savefig("figures/r3_vnia_nj.tiff", dpi=300, bbox_inches="tight")
    print("figures/r3_vnia_nj.pdf / .tiff")


if __name__ == "__main__":
    main(sys.argv[1])
