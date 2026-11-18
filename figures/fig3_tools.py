"""Fig. 3  cladeprint against an exact-k-mer baseline, Kover and Neptune (R4, 204-run subset).

A  share of each tool's markers (2,000 sampled per target) within 2 edits of a genome outside the target, among
   the unseen test genomes (census/out/bench/eval/summary.tsv, near2_test). Open circles count a marker when any
   one out-group genome is near it; filled circles when at least two are (census/out/bench/eval/nearmiss.tsv,
   ge2), which is the share left once the single tolerated out-group run is discounted. Neptune: VNII only.
B  wall time and peak memory on one 16-vCPU, 64 GB machine, each tool alone (census/out/bench/solo/*.time;
   Kover's memory from the container's anonymous-memory peak, logs/docker_peaks2.tsv, as /usr/bin/time
   measured only the docker client). Decimal GB.

Style as make_figures.py (SC Lab). Usage: micromamba run -n cladeprint python figures/fig3_tools.py
"""
import os
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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
TOOL = {"cp": ("cladeprint", "#2F4A3E", "o"), "cpx": ("exact k-mers", "#DFB166", "s"),
        "kover": ("Kover", "#C96A3A", "^"), "neptune": ("Neptune", "#9B2335", "D")}
B = "../census/out/bench"


def panel(ax, letter, dx=-0.2, dy=1.03):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")


def gnu_time(path):
    t = open(path).read()
    wall = re.search(r"Elapsed \(wall clock\).*: (\S+)", t).group(1).split(":")
    sec = sum(float(x) * 60 ** i for i, x in enumerate(reversed(wall)))
    kib = int(re.search(r"Maximum resident set size \(kbytes\): (\d+)", t).group(1))
    return sec / 60, kib * 1024 / 1e9


ev = pd.read_csv(f"{B}/eval/summary.tsv", sep="\t").dropna(subset=["near2_test"])
targets = list(dict.fromkeys(ev[ev.tool == "cp"].target))
short = {t: t.replace("C.neoformans/", "").replace("VNB/", "").replace("VNI/", "").replace("C.", "C. ")
         for t in targets}

fig, (a, b1, b2) = plt.subplots(1, 3, figsize=(7.0, 2.7), gridspec_kw={"width_ratios": [3.2, 1, 1]},
                                constrained_layout=True)
x = np.arange(len(targets))
# Neptune takes the exact-k-mer slot (those points sit near 1), so it does not cover cladeprint's circle
off = {"cpx": -0.2, "kover": 0.0, "cp": 0.2, "neptune": -0.2}
nm = pd.read_csv(f"{B}/eval/nearmiss.tsv", sep="\t")
nm = nm[nm.split == "test"].set_index("target")
for tool, (name, colour, m) in TOOL.items():
    g = ev[ev.tool == tool].set_index("target")
    xs = [i + off[tool] for i, t in enumerate(targets) if t in g.index]
    ys = [g.loc[t, "near2_test"] for t in targets if t in g.index]
    face = "none" if tool == "cp" else colour
    a.scatter(xs, ys, s=14 if tool != "neptune" else 18, marker=m, facecolors=face, edgecolors=colour
              if tool == "cp" else (SC_AXIS if tool == "neptune" else "none"), linewidths=0.8 if tool == "cp"
              else 0.4, color=None if tool == "cp" else colour,
              label={"cp": f"{name} (any one genome)", "neptune": f"{name} (VNII only)"}.get(tool, name), zorder=3)
    if tool == "cp":
        xs2 = [i + off[tool] for i, t in enumerate(targets) if t in nm.index]
        ys2 = [nm.loc[t, "ge2"] for t in targets if t in nm.index]
        a.scatter(xs2, ys2, s=14, marker=m, color=colour, edgecolors="none",
                  label="cladeprint (two or more)", zorder=4)
a.axvline(6.5, color="#E8E8E8", lw=0.75, zorder=0)
a.text(3, 1.04, "species", ha="center", fontsize=6.5, color=SC_ANNOT)
a.text(11, 1.04, "lineages and sublineages", ha="center", fontsize=6.5, color=SC_ANNOT)
a.set_xticks(x, [short[t] for t in targets], rotation=45, ha="right", fontsize=6.5)
for lab in a.get_xticklabels():
    if lab.get_text().startswith("C. "):
        lab.set_style("italic")
a.set(ylim=(0, 1.08), ylabel="Markers within 2 edits\nof an out-group test genome")
# legend in the empty lineage-side band, clear of every point
a.legend(loc="center", bbox_to_anchor=(0.66, 0.6), handletextpad=0.1)
panel(a, "A", dx=-0.09)

cp_min, cp_gb = gnu_time(f"{B}/solo/cladeprint_clean.time")
kv_min, _ = gnu_time(f"{B}/solo/kover_create.time")
kv_gb = int(open("../logs/docker_peaks2.tsv").read().split()[-1]) / 1e9
np_min, np_gb = gnu_time(f"{B}/solo/neptune_p2_C.neoformans_VNII.time")
rows = [("cladeprint\n(15 groups)", cp_min, cp_gb, "cp"), ("Kover\n(matrix only)", kv_min, kv_gb, "kover"),
        ("Neptune\n(VNII only)", np_min, np_gb, "neptune")]
y = np.arange(len(rows))[::-1]
for ax, k, xlabel in ((b1, 1, "Wall time (min)"), (b2, 2, "Peak memory (GB)")):
    ax.barh(y, [r[k] for r in rows], color=[TOOL[r[3]][1] for r in rows], alpha=0.8, height=0.6)
    for yi, r in zip(y, rows):
        ax.text(r[k] * 1.04, yi, f"{r[k]:.0f}" if r[k] >= 10 else f"{r[k]:.1f}", va="center", fontsize=6.5)
    ax.set_xlim(0, max(r[k] for r in rows) * 1.3)
    ax.set_xlabel(xlabel)
b1.set_yticks(y, [r[0] for r in rows], fontsize=6.5)
b2.set_yticks([])
b2.spines["left"].set_visible(False)
panel(b1, "B", dx=-0.95)
for ext, dpi in (("pdf", 600), ("png", 300), ("tiff", 300)):
    fig.savefig(f"fig3_tools.{ext}", dpi=dpi)
print(f"written: fig3_tools  (cladeprint {cp_min:.1f} min {cp_gb:.2f} GB; Kover {kv_min:.1f} min {kv_gb:.2f} GB;"
      f" Neptune {np_min:.1f} min {np_gb:.2f} GB)")
