"""Fig. S1  Parameter sensitivity (R6 sweep and the clone-cutoff rerun) against the census bundle.

A  marker count of each group relative to the default settings (f_in 0.90, max_out 1, clones at 5 SNPs/Mb).
B  5th-percentile leave-one-out sensitivity of each group, default and each alternative setting.
One setting changes at a time. Capped groups are compared on their stored hash samples (same rate in every run).

Inputs: census/out/{groups,loo}.tsv, census/out/sweep/{fin085,fin095,out0,out2}/, census/out/cut2/.
Style as make_figures.py (SC Lab). Usage: micromamba run -n cladeprint python figures/figS1_sweep.py
"""
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

os.chdir(os.path.dirname(os.path.abspath(__file__)))
SC_ANNOT, SC_AXIS, SC_GRID = "#4A4A4A", "#2C2C2C", "#E8E8E8"
mpl.rcParams.update({
    "font.family": ["Arial", "Liberation Sans", "DejaVu Sans"], "font.size": 7, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.75, "axes.edgecolor": SC_AXIS,
    "xtick.direction": "out", "ytick.direction": "out", "xtick.major.size": 3.5, "ytick.major.size": 3.5,
    "xtick.major.width": 0.75, "ytick.major.width": 0.75, "xtick.color": SC_AXIS, "ytick.color": SC_AXIS,
    "lines.linewidth": 1.5, "lines.markersize": 5, "legend.frameon": False, "figure.dpi": 300,
    "savefig.facecolor": "white", "savefig.bbox": "tight", "pdf.fonttype": 42, "ps.fonttype": 42,
})

RUNS = {"default": "../census/out",
        "f_in 0.85": "../census/out/sweep/fin085", "f_in 0.95": "../census/out/sweep/fin095",
        "max_out 0": "../census/out/sweep/out0", "max_out 2": "../census/out/sweep/out2",
        "clones at 2 SNPs/Mb": "../census/out/cut2"}
STYLE = {"f_in 0.85": ("#DFB166", "o"), "f_in 0.95": ("#9B2335", "o"),
         "max_out 0": ("#5E7A52", "s"), "max_out 2": ("#2F4A3E", "s"),
         "clones at 2 SNPs/Mb": ("#C96A3A", "^")}
ORDER = ["C.neoformans", "C.deneoformans", "C.gattii", "C.deuterogattii", "C.bacillisporus", "C.tetragattii",
         "C.decagattii", "C.neoformans/VNI", "C.neoformans/VNII", "C.neoformans/VNB",
         "C.neoformans/VNB/VNBI", "C.neoformans/VNB/VNBII", "C.neoformans/VNI/VNIa",
         "C.neoformans/VNI/VNIb", "C.neoformans/VNI/VNIc"]


def label(g):
    return g.rsplit("/", 1)[1] if "/" in g else "C. " + g.split(".", 1)[1]


def p05(v):
    v = sorted(v)
    return v[round((len(v) - 1) * 0.05)]


def panel(ax, letter, dx=-0.32, dy=1.02):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")


markers, sens = {}, {}
for k, d in RUNS.items():
    markers[k] = pd.read_csv(f"{d}/groups.tsv", sep="\t").set_index("group").markers
    sens[k] = pd.read_csv(f"{d}/loo.tsv", sep="\t").groupby("group").sensitivity.apply(p05)
markers, sens = pd.DataFrame(markers).loc[ORDER], pd.DataFrame(sens).loc[ORDER]
change = (markers.div(markers["default"], axis=0) - 1) * 100
change.round(2).join(sens.round(4), lsuffix=" markers %", rsuffix=" loo_p05").to_csv(
    "figS1_sweep.tsv", sep="\t")

y = range(len(ORDER))[::-1]
fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 3.3), sharey=True, constrained_layout=True)
for ax in (a, b):
    for yy, g in zip(y, ORDER):
        if g.count("/") == 1 and g.endswith("VNI"):
            ax.axhline(yy + 0.5, color=SC_GRID, lw=0.75, zorder=0)  # species | lineages
        if g.endswith("VNBI"):
            ax.axhline(yy + 0.5, color=SC_GRID, lw=0.75, zorder=0)  # lineages | sublineages
    ax.grid(axis="x", color=SC_GRID, lw=0.5, zorder=0)

# A: marker count change
a.axvline(0, color=SC_ANNOT, lw=0.75, zorder=1)
for k, (col, m) in STYLE.items():
    a.scatter(change[k], list(y), s=16, color=col, marker=m, linewidths=0, zorder=3, label=k)
a.set(yticks=list(y), yticklabels=[label(g) for g in ORDER], xlim=(-85, 85),
      xlabel="Markers relative to default (%)")
for t, g in zip(a.get_yticklabels(), ORDER):
    t.set_fontstyle("normal" if "/" in g else "italic")
panel(a, "A")

# B: 5th-percentile leave-one-out sensitivity
b.scatter(sens["default"], list(y), s=22, facecolors="none", edgecolors=SC_ANNOT, linewidths=0.9, zorder=4,
          label="default")
for k, (col, m) in STYLE.items():
    b.scatter(sens[k], list(y), s=16, color=col, marker=m, linewidths=0, zorder=3, label=k)
b.set(xlim=(0.4, 1.0), xlabel="Leave-one-out sensitivity, 5th percentile")
h, l = b.get_legend_handles_labels()
b.legend(h, l, loc="upper left", bbox_to_anchor=(0.02, 0.98), handletextpad=0.1)
panel(b, "B", dx=-0.04)

for ext in ("pdf", "png", "tiff"):
    fig.savefig(f"figS1_sweep.{ext}", dpi=600 if ext == "pdf" else 300)
print("written: figS1_sweep")
