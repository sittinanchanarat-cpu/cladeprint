# Figure S2 of the paper.
"""Fig. S2  Pure isolates left out of training, mostly at low coverage (census/ext/summary.tsv):
own-sublineage score against the share of the genome seen, coloured by the confidence of the call.

Style as make_figures.py (SC Lab). Usage: python figures/msystems/figS2_lowcov.py
"""
import os
import sys
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
os.makedirs(OUT, exist_ok=True)

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import glob as _glob
import matplotlib.font_manager as _fm
for _f in _glob.glob(os.path.join(os.environ.get("CLADEPRINT_FONTS", ""), "*.ttf")):
    _fm.fontManager.addfont(_f)
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))  # figures/; inputs are ../census/ext/
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
GROUP = {"VNI": ("C.neoformans/VNI", "#2F4A3E"), "VNII": ("C.neoformans/VNII", "#9B2335"),
         "VNB": ("C.neoformans/VNB", "#C96A3A"), "C. deneoformans": ("C.deneoformans", "#B8C4A2"),
         "C. gattii": ("C.gattii", "#DFB166")}
CONF_COLOR = {"high": "#2F4A3E", "medium": "#5E7A52", "low": "#DFB166", "none": "#D9D9D9"}


def panel(ax, letter, dx=-0.2, dy=1.03):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")


nn = pd.read_csv("../census/ext/nn_check.tsv", sep="\t", keep_default_na=False)
s = pd.read_csv("../census/ext/summary.tsv", sep="\t", keep_default_na=False)

NAME="FigureS2"
fig, c = plt.subplots(figsize=(3.4, 3.0), constrained_layout=True)
# C: untrained pure isolates, nearly all at low coverage
PARENT = {"VNBI": "VNB", "VNBII": "VNB", "VNIa": "VNI", "VNIb": "VNI", "VNIc": "VNI"}
p = s[s.kind == "pure"].copy()
p["own"] = [float(r[f"score:C.neoformans/{PARENT[r.label]}/{r.label}"]) for _, r in p.iterrows()]
p["conf"] = p.lineage_confidence.where(p.lineage != "NA", "none")
for conf in ["none", "low", "medium", "high"]:
    q = p[p.conf == conf]
    c.scatter(q.kmer_excess, q.own, s=14, color=CONF_COLOR[conf], edgecolors=SC_AXIS if conf == "none" else "none",
              linewidths=0.3, label={"none": "no lineage call"}.get(conf, conf), zorder=3)
c.plot([0, 1.05], [0, 1.05], color=SC_ANNOT, lw=0.75, ls="--", zorder=0)
c.axhline(0.3, color=SC_ANNOT, lw=0.75, ls=":")
c.set(xlim=(0, 1.05), ylim=(0, 1.05), xlabel="Share of genome k-mers seen", ylabel="Own-sublineage score")
c.legend(loc="upper left", handletextpad=0.1, title=f"{len(p)} isolates, 0 wrong", title_fontsize=6.5)

for ext, dpi in (("pdf", 600), ("png", 300), ("tiff", 300)):
    fig.savefig(os.path.join(OUT, f"{NAME}.{ext}"), dpi=dpi)
print("written", NAME)
