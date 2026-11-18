"""Fig. 1  Overview: the cladeprint workflow (A) and the census it was run on (B, C).

B  training runs per species after QC (census/out/qc.tsv, census/out/samples.tsv); excluded runs noted.
C  C. neoformans training runs by literature lineage label.

Style as make_figures.py (SC Lab). Usage: micromamba run -n cladeprint python figures/fig1_overview.py
"""
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

os.chdir(os.path.dirname(os.path.abspath(__file__)))
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
mpl.rcParams.update({"mathtext.fontset": "custom", "mathtext.it": "Arial:italic", "mathtext.rm": "Arial"})
DATA, METHOD, OUT = "#F5ECD7", "#B8C4A2", "#DFB166"  # box fills: data, marker discovery, outputs
LIN = {"VNI": "#2F4A3E", "VNB": "#C96A3A", "VNII": "#9B2335", "no lineage label": "#D9D9D9"}


def panel(ax, letter, dx, dy):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")


sm = pd.read_csv("../census/out/samples.tsv", sep="\t", dtype=str).fillna("").set_index("sample")
qc = pd.read_csv("../census/out/qc.tsv", sep="\t").set_index("sample")
d = qc.join(sm[["species", "lineage"]])
d["train"] = d.train.astype(str).str.lower() == "true"
tr = d[d.train]

fig = plt.figure(figsize=(7.0, 4.6), constrained_layout=True)
gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1], width_ratios=[1.35, 1])
a = fig.add_subplot(gs[0, :])
b, c = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])

# A: workflow, two rows read left to right, then the bundle feeds two outputs
a.set_xlim(0, 100)
a.set_ylim(0, 44)
a.axis("off")


def box(x, y, w, h, title, body, fill):
    a.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2", fc=fill, ec=SC_AXIS,
                               lw=0.6))
    a.text(x + w / 2, y + h - 2.2, title, ha="center", va="top", fontsize=7, fontweight="bold")
    a.text(x + w / 2, y + h - 6.4, body, ha="center", va="top", fontsize=6, linespacing=1.25)


def arrow(x0, y0, x1, y1, **kw):
    a.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=7, lw=0.8, color=SC_AXIS,
                                shrinkA=0, shrinkB=0, **kw))


w, h, top, low = 21, 17, 25, 2
row1 = [("Public reads", "3,736 runs, all seven\nspecies (NCBI SRA)", DATA),
        ("Counting", "31-mers per run (KMC),\nerror cutoff from the\nrun's own histogram", DATA),
        ("Labels", "species from the reads\n(331 NCBI names corrected);\nlineages from 3 studies", DATA),
        ("Clonal clusters", "5 SNPs/Mb; each cluster\nweighted as one isolate", METHOD)]
xs = [1, 26, 51, 76]
for x, (t, body, f) in zip(xs, row1):
    box(x, top, w, h, t, body, f)
for x in xs[:-1]:
    arrow(x + w + 0.8, top + h / 2, x + 25 - 0.8, top + h / 2)
row2 = [("Group markers", "in ≥ 90% of the group\nand of every child group;\nabsent from the out-group", METHOD),
        ("Near-miss filters", "no out-group k-mer within\n2 substitutions or\n2 edits (indels)", METHOD),
        ("Bundle", "15 groups at 3 levels:\nspecies, lineage,\nsublineage", METHOD)]
for x, (t, body, f) in zip(xs[:3], row2):
    box(x, low, w, h, t, body, f)
for x in xs[:2]:
    arrow(x + w + 0.8, low + h / 2, x + 25 - 0.8, low + h / 2)
ym = (top + low + h) / 2  # clusters -> markers: down, back to the left, down into the box
a.plot([xs[3] + w / 2, xs[3] + w / 2, xs[0] + w / 2], [top - 0.8, ym, ym], color=SC_AXIS, lw=0.8)
arrow(xs[0] + w / 2, ym, xs[0] + w / 2, low + h + 0.8)
a.add_patch(FancyBboxPatch((76, low + 9.2), w, 7.8, boxstyle="round,pad=0.4,rounding_size=1.2", fc=OUT,
                           ec=SC_AXIS, lw=0.6))
a.text(86.5, low + 15.4, "classify", ha="center", va="top", fontsize=7, fontweight="bold")
a.text(86.5, low + 11.8, "call, confidence, flags", ha="center", va="top", fontsize=6)
a.add_patch(FancyBboxPatch((76, low), w, 7.8, boxstyle="round,pad=0.4,rounding_size=1.2", fc=OUT, ec=SC_AXIS,
                           lw=0.6))
a.text(86.5, low + 6.2, "assays", ha="center", va="top", fontsize=7, fontweight="bold")
a.text(86.5, low + 2.6, "audit published; design new", ha="center", va="top", fontsize=6)
arrow(xs[2] + w + 0.8, low + h / 2, 75.2, low + 13)
arrow(xs[2] + w + 0.8, low + h / 2, 75.2, low + 4)
panel(a, "A", dx=-0.005, dy=0.98)

# B: training runs per species (log scale: C. neoformans dominates)
order = ["C.neoformans", "C.deuterogattii", "C.tetragattii", "C.gattii", "C.deneoformans", "C.bacillisporus",
         "C.decagattii"]
n = tr.species.value_counts().reindex(order)
y = range(len(order))[::-1]
b.barh(list(y), n.values, color="#5E7A52", alpha=0.8, height=0.65)
for yi, v in zip(y, n.values):
    b.text(v * 1.12, yi, f"{v:,}", va="center", fontsize=6.5)
b.set_xscale("log")
b.set_xlim(5, 12000)
b.set_yticks(list(y), [s.replace("C.", "C. ") for s in order], fontstyle="italic")
b.set_xlabel("Training runs")
ex = d[~d.train]
b.text(1.0, 0.02, f"{len(tr):,} training runs; {len(ex)} excluded\n({(ex.species == 'AD-hybrid').sum()} AD hybrids, "
       f"low coverage, mixtures, foreign)", transform=b.transAxes, ha="right", va="bottom", fontsize=6,
       color=SC_ANNOT)
panel(b, "B", dx=-0.33, dy=1.02)

# C: C. neoformans training runs by literature lineage label
neo = tr[tr.species == "C.neoformans"].lineage.replace("", "no lineage label").value_counts()
left = 0
for lin, colour in LIN.items():
    v = neo.get(lin, 0)
    c.barh(0, v, left=left, color=colour, alpha=0.8, height=0.5, edgecolor="white", linewidth=0.5)
    c.text(left + v / 2, 0.38 if lin != "VNII" else -0.38, f"{lin}\n{v:,}", ha="center",
           va="bottom" if lin != "VNII" else "top", fontsize=6.5)
    left += v
c.set_xlim(0, left)
c.set_ylim(-1, 1.6)
c.set_yticks([])
c.spines["left"].set_visible(False)
c.set_xlabel("$C.$ $neoformans$ training runs")
c.text(0.0, 0.97, "lineage labels: Rhodes 2017, Desjardins 2017,\nAshton 2019 (VNI sublineages: Ashton)",
       transform=c.transAxes, fontsize=6, color=SC_ANNOT, va="top")
panel(c, "C", dx=-0.04, dy=1.02)

for ext, dpi in (("pdf", 600), ("png", 300), ("tiff", 300)):
    fig.savefig(f"fig1_overview.{ext}", dpi=dpi)
print("written: fig1_overview")
