# Figure 4 of the paper.
"""Fig. 4  Classifying isolates the bundle was not trained on (census bundle).

A  lineage-unlabelled C. neoformans runs: SNP distance to the nearest literature-labelled run, by the
   lineage the markers called (census/ext/nn_check.tsv); every call agrees with that neighbour.
B  hybrid and diploid isolates: scores of the parent groups, with MIXED / EXCESS_KMERS flags
   (census/ext/summary.tsv).

Style as make_figures.py (SC Lab). Usage: python figures/msystems/fig4.py
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
    return ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")


nn = pd.read_csv("../census/ext/nn_check.tsv", sep="\t", keep_default_na=False)
s = pd.read_csv("../census/ext/summary.tsv", sep="\t", keep_default_na=False)
assert (nn.lineage == nn.nn_lineage).all(), "panel A says every call agrees"

NAME="Figure4"
fig = plt.figure(figsize=(7.0, 5.2), constrained_layout=True)
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.15])
a, b = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
# A: how far the unlabelled runs are from any labelled run, and that the calls agree anyway
bins = np.logspace(0, np.log10(3000), 36)
bottom = np.zeros(len(bins) - 1)
for lin in ["VNI", "VNB", "VNII"]:
    h, _ = np.histogram(nn[nn.lineage == lin].snps_per_mb.clip(lower=1), bins)
    a.bar(bins[:-1], h, np.diff(bins), bottom=bottom, align="edge", color=GROUP[lin][1], alpha=0.8,
          edgecolor="white", linewidth=0.2, label=f"{lin} ({(nn.lineage == lin).sum():,})")
    bottom += h
a.set_xscale("log")
a.set_yscale("symlog", linthresh=10)
a.set_yticks([0, 1, 10, 100, 1000], ["0", "1", "10", "100", "1,000"])
a.axvline(500, color=SC_ANNOT, lw=0.75, ls=":")
far = (nn.snps_per_mb > 500).sum()
a.text(560, 300, f"{far} runs\n> 500 SNPs/Mb", fontsize=6, color=SC_ANNOT, va="center")
a.set(xlabel="SNPs/Mb to nearest labelled run", ylabel="Unlabelled runs")
a.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, columnspacing=0.8, handlelength=0.9,
         handletextpad=0.3, title=f"lineage called; all {len(nn):,} agree with the neighbour", title_fontsize=6.5)
tA = panel(a, "A", dx=-0.21, dy=1.085)

# B: hybrids and diploids, one row each
order = ["AD hybrid", "A x gattii", "lineage hybrid", "diploid, one lineage"]
title = {"AD hybrid": "AD hybrids", "A x gattii": "A × C. gattii", "lineage hybrid": "VNI/VNII × VNB",
         "diploid, one lineage": "Diploids, one lineage"}
h = s[s.kind.isin(order)].copy()
h["o"] = h.kind.map(order.index)
h = h.sort_values(["o", "label", "isolate"]).reset_index(drop=True)
y = 0
ticks, names = [], []
for kind, g in h.groupby("o", sort=True):
    b.text(-0.02, y - 0.6, title[order[kind]], transform=b.get_yaxis_transform(), ha="right", va="center",
           fontsize=7, fontweight="bold")
    y += 0.6
    for r in g.itertuples():
        b.axhline(y, color="#E8E8E8", lw=0.5, zorder=0)
        for k, (name, (col, colour)) in enumerate(GROUP.items()):
            v = float(h.loc[r.Index, f"score:{col}"])
            if v >= 0.05:
                b.scatter(v, y + (k - 2) * 0.1, s=16, color=colour, edgecolors=SC_AXIS, linewidths=0.3, zorder=3)
        mark = ("M" if r.mixed else "") + ("E" if r.excess else "")
        b.text(1.03, y, mark, transform=b.get_yaxis_transform(), fontsize=6, va="center", color=SC_ANNOT)
        ticks.append(y)
        names.append(f"{r.isolate}  ({r.label.replace('VNIV', 'D').replace('VG', 'gattii')})")
        y += 1
    y += 0.4
b.set_yticks(ticks, names, fontsize=6)
b.set_ylim(y - 0.6, -0.8)
b.set_xlim(0, 1)
b.axvline(0.3, color=SC_ANNOT, lw=0.75, ls=":")
b.set_xlabel("Score (parent groups)")
for name, (_, colour) in GROUP.items():
    b.scatter([], [], s=16, color=colour, edgecolors=SC_AXIS, linewidths=0.3, label=name)
b.legend(loc="lower left", bbox_to_anchor=(-0.05, 1.0), ncol=3, columnspacing=0.6, handletextpad=0.1)
b.text(1.03, 1.0, "M  MIXED\nE  EXCESS_KMERS", transform=b.transAxes, fontsize=6, va="bottom", color=SC_ANNOT)
tB = panel(b, "B", dx=-0.62, dy=1.075)

for ext, dpi in (("pdf", 600), ("png", 300), ("tiff", 300)):
    fig.savefig(os.path.join(OUT, f"{NAME}.{ext}"), dpi=dpi)
print("written", NAME)
