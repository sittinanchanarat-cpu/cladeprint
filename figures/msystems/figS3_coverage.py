# Figure S3 of the paper.
"""Fig. S3  Coverage calibration (census/calib/summary.tsv from census/calib/summarize.py).

A  score of the true group vs share of the genome seen (kmer_excess), all three levels.
B  outcome by base coverage, per level: high / medium / low confidence / no call (no wrong calls).
C  genome k-mer peak vs base coverage, species level, coloured by confidence.

Style as make_figures.py (SC Lab). Usage: python figures/msystems/figS3_coverage.py
"""
import os

import matplotlib as mpl
import glob as _glob
import matplotlib.font_manager as _fm
for _f in _glob.glob(os.path.join(os.environ.get("CLADEPRINT_FONTS", ""), "*.ttf")):
    _fm.fontManager.addfont(_f)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))  # repo root (figures/msystems/ -> repo)
os.chdir(HERE)
os.makedirs("out", exist_ok=True)
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
LEVEL_COLOR = {"species": "#9B2335", "lineage": "#C96A3A", "sublineage": "#5E7A52"}
CONF_COLOR = {"high": "#2F4A3E", "medium": "#5E7A52", "low": "#DFB166", "none": "#D9D9D9"}
CONF_LABEL = {"high": "high", "medium": "medium", "low": "low", "none": "no call"}


def panel(ax, letter, dx=-0.2, dy=1.04):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")


s = pd.read_csv(os.path.join(ROOT, "census/calib/summary.tsv"), sep="\t")
assert (s.outcome != "wrong").all(), "caption says no wrong calls"

fig, (b, a, c) = plt.subplots(1, 3, figsize=(7.0, 2.5), gridspec_kw={"width_ratios": [1.25, 1, 1]},
                              constrained_layout=True)

# A: score tracks the share of the genome seen
for level in ["species", "lineage", "sublineage"]:
    g = s[(s.level == level) & (s.outcome == "right")]
    a.scatter(g.kmer_excess, g.score, s=7, color=LEVEL_COLOR[level], alpha=0.8, linewidths=0, label=level)
a.plot([0, 1.05], [0, 1.05], color=SC_ANNOT, lw=0.75, ls="--", zorder=0)
a.axhline(0.3, color=SC_ANNOT, lw=0.75, ls=":")
a.text(1.05, 0.31, "call threshold", ha="right", va="bottom", fontsize=6, color=SC_ANNOT)
a.set(xlim=(0, 1.07), ylim=(0, 1.05), xlabel="Share of genome k-mers seen", ylabel="Score of the true group")
a.legend(loc="upper left", handletextpad=0.1)
panel(a, "B")

# B: outcome by base coverage, three bars (levels) per bin
bins = [1, 2, 4, 8, 16, 32, 70]
names = ["1–2", "2–4", "4–8", "8–16", "16–32", "32–66"]
s["bin"] = pd.cut(s.coverage, bins, labels=names, right=False)
w = 0.27
for k, level in enumerate(["species", "lineage", "sublineage"]):
    g = s[s.level == level]
    share = pd.crosstab(g.bin, g.confidence, normalize="index").reindex(names).fillna(0)
    bottom = np.zeros(len(names))
    x = np.arange(len(names)) + (k - 1) * (w + 0.02)
    for conf in ["high", "medium", "low", "none"]:
        v = share[conf].values if conf in share else np.zeros(len(names))
        b.bar(x, v, w, bottom=bottom, color=CONF_COLOR[conf], alpha=0.8, edgecolor="white", linewidth=0.3,
              label=CONF_LABEL[conf] if k == 0 else None)
        bottom += v
    b.text(x[0], 0.5, level, ha="center", va="center", rotation=90, fontsize=6, color=SC_AXIS)
b.set_xticks(np.arange(len(names)), names)
b.set(ylim=(0, 1.02), xlabel="Base coverage (×)", ylabel="Share of run × depth")
b.set_yticks([0, 0.25, 0.5, 0.75, 1])
b.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), handlelength=0.9, title="confidence",
         title_fontsize=7)
panel(b, "A", dx=-0.17)

# C: the genome peak, not base coverage, separates high from lower confidence
g = s[s.level == "species"]
rng = np.random.default_rng(1)
for conf in ["none", "low", "medium", "high"]:
    h = g[g.confidence == conf]
    c.scatter(h.coverage, h.peak + rng.uniform(-0.25, 0.25, len(h)), s=8, color=CONF_COLOR[conf],
              edgecolors=SC_AXIS if conf == "none" else "none", linewidths=0.3, alpha=0.9, zorder=3)
c.axhline(4.5, color=SC_ANNOT, lw=0.75, ls=":")
c.text(1.1, 4.9, "peak ≥ 5: all high", fontsize=6, color=SC_ANNOT, va="bottom")
c.set_xscale("log")
c.set_xticks([1, 2, 4, 8, 16, 32, 64], ["1", "2", "4", "8", "16", "32", "64"])
c.minorticks_off()
c.set_yscale("symlog", linthresh=8)
c.set_yticks([0, 2, 4, 6, 8, 16, 32, 64], ["0", "2", "4", "6", "8", "16", "32", "64"])
c.set(ylim=(-0.6, 70), xlabel="Base coverage (×)", ylabel="Genome k-mer peak (species)")
panel(c, "C", dx=-0.2)

for ext, dpi in (("pdf", 600), ("png", 300), ("tiff", 300)):
    fig.savefig(f"out/FigureS3.{ext}", dpi=dpi)
print("written: figures/msystems/out/FigureS3")
