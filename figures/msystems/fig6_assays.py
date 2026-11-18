# Figure 6 of the paper.
"""Fig. 6  Published assays on public genomes, and new assays for groups that have none.

A  in-silico PCR of published assays on 274 NCBI assemblies (assay_panel.py, main-text rows).
B  product sizes reconstructed from the assemblies against the sizes the papers give
   (census/asm/r5_pcr_products.tsv, labels/assays/r5_pcr_reactions.tsv). Only products within +-10% (>= 20 bp)
   of a published size are counted, so the panel shows agreement within that window.
C  member vs non-member rate of the best designed allele-specific pair per group (3'-anchored on markers;
   census/out/design3p/pairs_screened.tsv), with the Kelley et al. 2014 MAMA assays on the same genomes
   (census/out/r5_kelley_mama.tsv, from census/r5_kelley_mama.py). Both sides are rates over the 3,426
   training runs grouped by second-round label: primer_sites screened 3,624 genomes, but
   design3p_screen_summary.py and r5_kelley_mama.py both restrict the rates to trained runs.
   Non-member rates below `floor` are plotted at `floor`; the stored values are unfloored.

Style as make_figures.py (SC Lab). Usage: python figures/msystems/fig6_assays.py
"""
import importlib.util
import os

import matplotlib as mpl
import glob as _glob
import matplotlib.font_manager as _fm
for _f in _glob.glob(os.path.join(os.environ.get("CLADEPRINT_FONTS", ""), "*.ttf")):
    _fm.fontManager.addfont(_f)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
os.makedirs("figures/msystems/out", exist_ok=True)
spec = importlib.util.spec_from_file_location("r5a", "figures/msystems/assay_panel.py")
r5a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(r5a)  # sets its own rcParams; ours follow
SC_ANNOT, SC_AXIS = "#4A4A4A", "#2C2C2C"
mpl.rcParams.update({
    "font.family": ["Arial", "Liberation Sans", "DejaVu Sans"], "font.size": 7, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "xtick.direction": "out", "ytick.direction": "out", "xtick.major.size": 3.5, "ytick.major.size": 3.5,
    "xtick.major.width": 0.75, "ytick.major.width": 0.75, "xtick.color": SC_AXIS, "ytick.color": SC_AXIS,
    "legend.frameon": False, "figure.dpi": 300, "savefig.facecolor": "white", "savefig.bbox": "tight",
    "pdf.fonttype": 42, "ps.fonttype": 42,
})
_kelley = pd.read_csv("census/out/r5_kelley_mama.tsv", sep="\t")
KELLEY = {r.target_species: (r.molecular_type, r.member, r.nonmember) for r in _kelley.itertuples()}


def panel(ax, letter, dx, dy=1.02):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")


def spines(ax):
    ax.spines[["top", "right"]].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(0.75)
        ax.spines[s].set_color(SC_AXIS)


fig = plt.figure(figsize=(7.0, 7.4), constrained_layout=True)
gs = fig.add_gridspec(2, 2, height_ratios=[1.55, 1])
a = fig.add_subplot(gs[0, :])
b, c = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])

# A: the published-assay heatmap
r5a.main(ax=a)
a.tick_params(length=0)
panel(a, "A", dx=-0.42, dy=1.05)

# B: reconstructed vs published product sizes
prod = pd.read_csv("census/asm/r5_pcr_products.tsv", sep="\t")
prod = prod[prod.expected & prod.probe_ok & (prod.rule == "lenient")]
rx = pd.read_csv("labels/assays/r5_pcr_reactions.tsv", sep="\t").dropna(subset=["expected_bp"])
pts = []
for r in rx.itertuples():
    exp = [int(x) for x in str(r.expected_bp).split(",")]
    p = prod[prod.reaction == r.reaction]
    for e in exp:
        s = p[(p["size"] - e).abs() <= max(20, 0.1 * e)]["size"]
        if len(s) >= 3:
            pts.append((e, float(np.median(s)), len(s), r.reaction))
pts = pd.DataFrame(pts, columns=["published", "predicted", "n", "reaction"])
b.plot([50, 3500], [50, 3500], color=SC_ANNOT, lw=0.75, ls="--", zorder=0)
b.scatter(pts.published, pts.predicted, s=8 + 2 * np.sqrt(pts.n), color="#C96A3A", alpha=0.8,
          edgecolors=SC_AXIS, linewidths=0.3, zorder=3)
b.set(xscale="log", yscale="log", xlim=(50, 3500), ylim=(50, 3500), xlabel="Published product (bp)",
      ylabel="Reconstructed product (bp, median)")
dev = (pts.predicted / pts.published - 1).abs()
b.text(0.04, 0.96, f"{len(pts)} products in {pts.reaction.nunique()} reactions\nmedian deviation "
       f"{100 * dev.median():.1f}%", transform=b.transAxes, va="top", fontsize=6.5)
spines(b)
panel(b, "B", dx=-0.2)

# C: designed pairs vs Kelley MAMA
d = pd.read_csv("census/out/design3p/pairs_screened.tsv", sep="\t")
best = d.sort_values("margin", ascending=False).groupby("group").head(1).set_index("group")
LEVEL = lambda g: "species" if "/" not in g else ("lineage / sublineage" if g.count("/") <= 2 else "VNIa sub-clade")
COL = {"species": "#2F4A3E", "lineage / sublineage": "#C96A3A", "VNIa sub-clade": "#B8C4A2"}
floor = 3e-4
for lv, colour in COL.items():
    g = best[best.index.map(LEVEL) == lv]
    c.scatter(g.nonmember.clip(lower=floor), g.member, s=16, color=colour, edgecolors=SC_AXIS, linewidths=0.3,
              label=f"designed: {lv.replace(' / ', '/')}",
              zorder=3.5 if lv == "lineage / sublineage" else 3)  # VNII sits between two sub-clade points
for grp, (vg, mem, non) in KELLEY.items():
    c.scatter(max(non, floor), mem, s=22, marker="D", facecolors="none", edgecolors="#9B2335", linewidths=0.9,
              zorder=4, label="Kelley 2014 MAMA" if vg == "VGI" else None)
# VGI, VGII and VGIV all sit at the floor with near-equal member shares, on top of the species points:
# label them in a column in the empty band below, joined by leader lines kept in the same order.
floor_vg = sorted((v for v in KELLEY.values() if v[2] < floor), key=lambda v: -v[1])
for k, (vg, mem, non) in enumerate(floor_vg):
    c.annotate(vg, (floor, mem), xytext=(9e-4, 0.9885 - 0.0045 * k), textcoords="data", fontsize=5.5,
               color="#9B2335", va="center",
               arrowprops=dict(arrowstyle="-", lw=0.4, color="#9B2335", shrinkA=1, shrinkB=3))
for vg, mem, non in KELLEY.values():
    if non >= floor:
        c.annotate(vg, (non, mem), xytext=(-6, -1), textcoords="offset points", fontsize=5.5,
                   color="#9B2335", ha="right", va="center")
# Every lineage/sublineage pair is labelled, each label placed in empty space (data coordinates) with a
# leader line, so no label sits between two points; VNIa and VNIc are clear of neighbours and label directly.
YTOP = 1.012
for grp, (tx, ty, ha) in {"C.neoformans/VNB/VNBI": (5.2e-3, 0.9855, "center"),
                          "C.neoformans/VNB/VNBII": (1.25e-2, 0.9695, "left"),
                          "C.neoformans/VNB": (2.8e-2, 0.9785, "left"),
                          "C.neoformans/VNII": (1.44e-2, 1.0065, "center"),
                          "C.neoformans/VNI/VNIb": (1.24e-1, 1.0065, "center"),
                          "C.neoformans/VNI": (6.2e-1, 0.983, "center")}.items():
    r = best.loc[grp]
    c.annotate(grp.split("/")[-1], (max(r.nonmember, floor), r.member), xytext=(tx, ty), textcoords="data",
               fontsize=5.5, ha=ha, va="center",
               arrowprops=dict(arrowstyle="-", lw=0.4, color=SC_ANNOT, shrinkA=1, shrinkB=3))
r = best.loc["C.neoformans/VNI/VNIa"]
c.annotate("VNIa", (r.nonmember, r.member), xytext=(5, -3), textcoords="offset points", fontsize=5.5)
c.axvspan(floor * 0.7, 0.01, ymin=(0.95 - 0.9) / (YTOP - 0.9), color="#E8E8E8", zorder=0, lw=0)
c.text(floor * 0.85, 0.9515, "member ≥ 0.95, non-member ≤ 0.01", ha="left", va="bottom", fontsize=5.5,
       color=SC_ANNOT)
c.annotate("VNIc", (best.loc["C.neoformans/VNI/VNIc", "nonmember"], best.loc["C.neoformans/VNI/VNIc", "member"]),
           xytext=(-5, 0), textcoords="offset points", fontsize=5.5, ha="right", va="center")
c.set(xscale="log", xlim=(floor * 0.7, 1), ylim=(0.9, YTOP), xlabel="Non-member genomes amplified",
      ylabel="Member genomes amplified")
c.set_xticks([1e-3, 1e-2, 1e-1, 1], ["≤0.001", "0.01", "0.1", "1"])
c.legend(loc="lower right", fontsize=5.5, handletextpad=0.1, borderaxespad=0.2)
spines(c)
panel(c, "C", dx=-0.2)

for ext, dpi in (("pdf", 600), ("png", 300), ("tiff", 300)):
    fig.savefig(f"figures/msystems/out/Figure6.{ext}", dpi=dpi)
print("written: figures/msystems/out/Figure6;", len(pts), "products; median dev", round(100 * dev.median(), 2), "%")
