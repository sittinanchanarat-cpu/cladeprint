"""Paper figures (SC Lab style). Reads the result tables; writes figures/fig*.{pdf,png,tiff}.

Fig. 2  Marker discovery: C. neoformans species markers carried per lineage, flat vs hierarchy-aware
        prevalence (A, census); attrition from exact group-specific 31-mers through the Hamming and edit
        near-miss filters (B, k-mers; C, independent loci on H99; pilot).

Usage: micromamba run -n cladeprint python figures/make_figures.py
"""
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

os.chdir(os.path.dirname(os.path.abspath(__file__)))
R = ".."

# ── SC style (sc-figures skill); Liberation Sans is metric-compatible with Arial ──────────────
SC_ANNOT, SC_AXIS = "#4A4A4A", "#2C2C2C"
SC_WARM_CMAP = mpl.colors.LinearSegmentedColormap.from_list("sc_warm", ["#F5ECD7", "#DFB166", "#C96A3A", "#9B2335"])
mpl.rcParams.update({
    "font.family": ["Arial", "Liberation Sans", "DejaVu Sans"], "font.size": 7, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.75, "axes.edgecolor": SC_AXIS,
    "xtick.direction": "out", "ytick.direction": "out", "xtick.major.size": 3.5, "ytick.major.size": 3.5,
    "xtick.major.width": 0.75, "ytick.major.width": 0.75, "xtick.color": SC_AXIS, "ytick.color": SC_AXIS,
    "lines.linewidth": 1.5, "lines.markersize": 5, "legend.frameon": False, "figure.dpi": 300,
    "savefig.facecolor": "white", "savefig.bbox": "tight", "pdf.fonttype": 42, "ps.fonttype": 42,
})
TARGETS = ["VNII", "VNBI", "VNBII", "VNI", "VNIa", "VNIb", "VNIc"]
COLOR = {"VNII": "#9B2335", "VNBI": "#C96A3A", "VNBII": "#DFB166", "VNI": "#2F4A3E",
         "VNIa": "#5E7A52", "VNIb": "#B8C4A2", "VNIc": SC_ANNOT}
SUB = {"VNIa": "VNI", "VNIb": "VNI", "VNIc": "VNI"}


def panel(ax, letter, dx=-0.16, dy=1.04):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")


def save(fig, name):
    for ext, dpi in (("pdf", 600), ("png", 300), ("tiff", 300)):
        fig.savefig(f"{name}.{ext}", dpi=dpi)
    plt.close(fig)


# ── Fig. 2: prevalence in every child group, and marker attrition ────────────────────────────
nm = pd.read_csv(f"{R}/kmers/near_miss.tsv", sep="\t").set_index("target")
ed = pd.read_csv(f"{R}/kmers/edit_miss.tsv", sep="\t").set_index("target")
pl = pd.read_csv(f"{R}/kmers/place_markers.tsv", sep="\t")
loci = pl.pivot(index="target", columns="set", values="loci")
lab = pd.read_csv(f"{R}/census/out/samples.tsv", sep="\t", dtype=str).fillna("").set_index("sample")


def species_scores(f):
    t = pd.read_csv(f, sep="\t").set_index("sample")
    t = t[t.train.astype(str).str.lower() == "true"]
    t = t[t.index.map(lab.species) == "C.neoformans"]
    return t.assign(lin=t.index.map(lab.lineage))


flat = species_scores(f"{R}/census/out/scores_flat_v1.tsv")
bal = species_scores(f"{R}/census/out/scores.tsv")

fig, (p, a, b) = plt.subplots(1, 3, figsize=(7.0, 2.9), gridspec_kw={"width_ratios": [0.85, 1.25, 1]},
                              constrained_layout=True)
# A: C. neoformans species markers carried by each lineage's training runs, flat vs hierarchy-aware
rng = np.random.default_rng(1)
LIN = ["VNI", "VNB", "VNII"]
for i, lin in enumerate(LIN):
    for k, (t, name) in enumerate([(flat, "flat"), (bal, "every child")]):
        v = t[t.lin == lin]["C.neoformans"].values
        x = i + (k - 0.5) * 0.38
        sub = v if len(v) <= 200 else rng.choice(v, 200, replace=False)  # thin VNI for legibility
        p.scatter(x + rng.uniform(-0.1, 0.1, len(sub)), sub, s=3, color=COLOR[lin] if lin != "VNB" else "#C96A3A",
                  alpha=0.35 if k == 0 else 0.8, linewidths=0, zorder=2)
        p.plot([x - 0.14, x + 0.14], [np.median(v)] * 2, color=SC_AXIS, lw=1.2, zorder=3)
        if lin == "VNII":
            p.text(x, np.median(v) - 0.04, f"{np.median(v):.2f}", ha="center", va="top", fontsize=6)
p.set_xticks(range(3), LIN)
for i in range(3):
    for k, lab_ in enumerate(["flat", "child"]):
        p.text(i + (k - 0.5) * 0.38, 0.33, lab_, ha="center", fontsize=5.5, color=SC_ANNOT)
p.set_ylim(0.3, 1.02)
p.set_ylabel("C. neoformans species markers present")
panel(p, "A", dx=-0.3)

# B: attrition through the near-miss filters
stages = ["Exact\n(absent from out-group)", "Hamming ≤ 2\nexcluded", "Edit ≤ 2\nexcluded"]
for t in TARGETS:
    y = [nm.loc[t, "strict"], nm.loc[t, "survive"], ed.loc[t, "survive"]]
    ls = "--" if t in SUB else "-"
    a.plot(range(3), y, ls=ls, marker="o", color=COLOR[t], ms=4, lw=1.2)
    nudge = {"VNIb": 1.12, "VNIa": 1.0, "VNIc": 0.86}.get(t, 1.0)  # keep the three small labels apart
    a.text(2.08, y[-1] * nudge, f"{t}  {y[-1]:,}", va="center", fontsize=6.5,
           color=COLOR[t] if t != "VNIb" else "#7D8B67")
a.set_yscale("log")
a.set_xticks(range(3), stages, fontsize=6.5)
a.set_xlim(-0.15, 2.9)
a.set_ylabel("Lineage-specific 31-mers")
a.spines["bottom"].set_bounds(0, 2)
panel(a, "B", dx=-0.17)

# C: independent loci left
y = np.arange(len(TARGETS))[::-1]
b.barh(y, loci.loc[TARGETS, "ed2"], color=[COLOR[t] for t in TARGETS], alpha=0.8, height=0.65,
       label="after edit filter")
b.scatter(loci.loc[TARGETS, "nm2"], y, facecolors="white", edgecolors=SC_AXIS, s=14, lw=0.75, zorder=3,
          label="after Hamming filter")
for yi, t in zip(y, TARGETS):
    v = loci.loc[t, "ed2"]
    inside = v >= 40
    b.text(v * (0.85 if inside else 1.15), yi, f"{int(v):,}", ha="right" if inside else "left", va="center",
           fontsize=6, color="white" if inside and t in ("VNII", "VNI") else SC_AXIS)
b.set_xscale("log")
b.set_yticks(y, TARGETS)
b.set_xlabel("Independent loci on H99")
b.legend(loc="lower right", handletextpad=0.4)
panel(b, "C", dx=-0.24)
save(fig, "fig2_markers")

print("written: fig2_markers")
