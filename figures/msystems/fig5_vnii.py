# Figure 5 of the paper.
"""Fig. 5  VNII carries a C. gattii-type gene cluster that the rest of C. neoformans lost.

A  the locus in C. gattii WM276 (chr F), VNII (from the long-read genome GCA_022832995.1 and five read sets;
   vnii/) and H99 (VNI; VNB is the same), aligned on the start of HOC1 / CGB_F4610W. Gene coordinates:
   vnii/wm276.gff, H99 from miniprot (vnii/genes/homologs.tsv); the KDZ insert length from the long-read genome.
B  the nine WM276 cluster genes in each species (miniprot % amino-acid identity; blank = absent or < 50% of the
   protein aligned; * frameshift or stop). VNII column: reads of five isolates at single-copy depth on every
   WM276-type gene and none on the H99-type ones (vnii/locus/locus_genes.txt).
C  maximum-likelihood tree of the 23-kb block (vnii/tree/block.contree; IQ-TREE, UFBoot), rooted on
   C. deneoformans.

Style as make_figures.py (SC Lab). Usage: python figures/msystems/fig5_vnii.py
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
from Bio import Phylo
from matplotlib.patches import Polygon, Rectangle

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
os.makedirs("figures/msystems/out", exist_ok=True)
SC_ANNOT, SC_AXIS = "#4A4A4A", "#2C2C2C"
mpl.rcParams.update({
    "font.family": ["Arial", "Liberation Sans", "DejaVu Sans"], "font.size": 7, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.75, "axes.edgecolor": SC_AXIS,
    "xtick.direction": "out", "ytick.direction": "out", "xtick.major.size": 3.5, "ytick.major.size": 3.5,
    "xtick.major.width": 0.75, "ytick.major.width": 0.75, "xtick.color": SC_AXIS, "ytick.color": SC_AXIS,
    "lines.linewidth": 1.5, "lines.markersize": 5, "legend.frameon": False, "figure.dpi": 300,
    "savefig.facecolor": "white", "savefig.bbox": "tight", "pdf.fonttype": 42, "ps.fonttype": 42,
    "mathtext.fontset": "custom", "mathtext.it": "Arial:italic", "mathtext.rm": "Arial",
})
GATTII, NEO, LOST, FLANK, KDZ = "#C96A3A", "#2F4A3E", "#9B2335", "#B0B0B0", "#DFB166"
BLOCK = ["CGB_F4610W", "CGB_F4620C", "CGB_F4630C", "CGB_F4640C", "CGB_F4650W", "CGB_F4660C", "CGB_F4670W",
         "CGB_F4680W", "CGB_F4690W"]
SHORT = {"CGB_F4610W": "HOC1", "CGB_F4620C": "4620", "CGB_F4630C": "4630", "CGB_F4640C": "4640",
         "CGB_F4650W": "4650", "CGB_F4660C": "4660", "CGB_F4670W": "4670", "CGB_F4680W": "4680",
         "CGB_F4690W": "RGD1"}
LOST_IN_NEO = {"CGB_F4630C", "CGB_F4640C", "CGB_F4660C", "CGB_F4670W", "CGB_F4680W"}


def panel(ax, letter, dx, dy=1.02):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")


def gene(ax, x0, x1, strand, y, colour, h=0.32, ls="-", fill=True):
    head = min(0.45, (x1 - x0) * 0.4)
    if strand == "+":
        pts = [(x0, y - h / 2), (x1 - head, y - h / 2), (x1, y), (x1 - head, y + h / 2), (x0, y + h / 2)]
    else:
        pts = [(x1, y - h / 2), (x0 + head, y - h / 2), (x0, y), (x0 + head, y + h / 2), (x1, y + h / 2)]
    ax.add_patch(Polygon(pts, closed=True, fc=colour if fill else "white", ec=SC_AXIS, lw=0.5, ls=ls, zorder=3))


gff = pd.read_csv("vnii/wm276.gff", sep="\t", comment="#", header=None)
gff = gff[(gff[2] == "gene") & (gff[0] == gff[0][gff[8].str.contains("CGB_F4610W")].iloc[0])]
gff["name"] = gff[8].str.extract(r"Name=([^;]+)")
wm = gff.set_index("name").loc[["CGB_F4600C"] + BLOCK + ["CGB_F4700C"], [3, 4, 6]]
wm.columns = ["s", "e", "strand"]
wm["s"], wm["e"] = wm.s.astype(float) / 1000, wm.e.astype(float) / 1000
ho = pd.read_csv("vnii/genes/homologs.tsv", sep="\t")
h99 = ho[(ho.species == "neoformans") & (ho.contig == "NC_026751.1")].set_index("gene")

fig = plt.figure(figsize=(7.0, 6.6), constrained_layout=True)
gs = fig.add_gridspec(3, 2, height_ratios=[1, 1.35, 0.78], width_ratios=[1.05, 1])
a = fig.add_subplot(gs[0, :])
b, c = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
d = fig.add_subplot(gs[2, :])

# A: three tracks, x in kb from the start of HOC1 / F4610W
o_wm, o_h99 = wm.loc["CGB_F4610W", "s"], h99.loc["CGB_F4610W", "start"] / 1000
rows = {"$C.$ $gattii$ WM276 (VGI)": 2.0, "VNII": 1.0, "H99 (VNI; VNB alike)": 0.0}
for label, y in rows.items():
    a.text(-24.5, y, label, ha="left", va="center", fontsize=7)
# WM276
a.plot([-6.5, 24], [2, 2], color=SC_AXIS, lw=0.6, zorder=1)
for g, r in wm.iterrows():
    gene(a, r.s - o_wm, r.e - o_wm, r.strand, 2.0, GATTII if g in BLOCK else FLANK)
# VNII: H99-type flank, KDZ insert (10.8 kb, one TIR), gattii-type block from 1.8 kb upstream of HOC1
blk0, blk1 = 923.314 - o_wm, 947.3 - o_wm
ins0 = blk0 - 10.834
a.plot([ins0 - 4.5, blk1 + 2.5], [1, 1], color=SC_AXIS, lw=0.6, zorder=1)
a.add_patch(Rectangle((ins0, 0.86), 10.834, 0.28, fc=KDZ, ec=SC_AXIS, lw=0.5, hatch="////", zorder=2))
a.text(ins0 + 5.4, 1.25, "KDZ transposon (10.8 kb)", ha="center", va="bottom", fontsize=6)
gene(a, ins0 - 4.2, ins0 - 2.3, "-", 1.0, FLANK)
for g in BLOCK:
    r = wm.loc[g]
    gene(a, r.s - o_wm, r.e - o_wm, r.strand, 1.0, GATTII)
gene(a, blk1 + 0.3, blk1 + 1.1, "-", 1.0, FLANK)
a.add_patch(Rectangle((blk0, 0.72), blk1 - blk0, 0.56, fc="none", ec=GATTII, lw=0.8, ls="--", zorder=2))
# H99
a.plot([-6.5, 9], [0, 0], color=SC_AXIS, lw=0.6, zorder=1)
for g, r in h99.iterrows():
    x0, x1 = r.start / 1000 - o_h99, r.end / 1000 - o_h99
    strand = wm.loc[g, "strand"] if g in wm.index else "+"
    if g == "CGB_F4650W":
        gene(a, x0, x1, strand, 0.0, LOST, ls="--", fill=False)
    else:
        gene(a, x0, x1, strand, 0.0, NEO if g in BLOCK else FLANK)
for k, g in enumerate(BLOCK):
    r = wm.loc[g]
    a.text((r.s + r.e) / 2 - o_wm, 2.27 + 0.22 * (k == 1), SHORT[g], ha="center", va="bottom", fontsize=5.5)
a.text(h99.loc["CGB_F4610W", "start"] / 1000 - o_h99 + 0.8, -0.32, "HOC1", ha="center", va="top", fontsize=5.5)
a.text(h99.loc["CGB_F4690W", "start"] / 1000 - o_h99 + 1.5, -0.32, "RGD1", ha="center", va="top", fontsize=5.5)
a.text(h99.loc["CGB_F4650W", "start"] / 1000 - o_h99 + 0.6, 0.3, "relic", ha="center", va="bottom", fontsize=5.5,
       color=LOST)
a.text(9.6, 0, "4630, 4640, 4660, 4670, 4680 absent", va="center", fontsize=6, color=LOST)
a.plot([0, 5], [-0.75, -0.75], color=SC_AXIS, lw=0.8)
a.text(2.5, -0.82, "5 kb", ha="center", va="top", fontsize=6)
a.set_xlim(-25, 26)
a.set_ylim(-1.05, 2.9)
a.axis("off")
panel(a, "A", dx=0.0, dy=0.96)

# B: gene presence and identity across species
cols = [("H99\n(VNI, VNB)", "neoformans"), ("VNII", None), ("C. deneo-\nformans", "deneoformans"),
        ("C. gattii", "gattii"), ("C. deutero-\ngattii", "deuterogattii"), ("C. bacilli-\nsporus", "bacillisporus"),
        ("C. tetra-\ngattii", "tetragattii"), ("C. deca-\ngattii", "decagattii")]
cmap = mpl.colors.LinearSegmentedColormap.from_list("sc_warm", ["#F5ECD7", "#DFB166", "#C96A3A", "#9B2335"])
norm = mpl.colors.Normalize(70, 100)
for j, (name, sp) in enumerate(cols):
    for i, g in enumerate(BLOCK):
        if sp is None:  # VNII: gattii-type copy present at single-copy depth
            b.add_patch(Rectangle((j - 0.45, i - 0.45), 0.9, 0.9, fc=GATTII, ec="white", lw=0.5))
            b.text(j, i, "+", ha="center", va="center", fontsize=7, color="white", fontweight="bold")
            continue
        r = ho[(ho.species == sp) & (ho.gene == g)]
        ok = len(r) and r.coverage.iloc[0] >= 0.5 and r.identity.iloc[0] >= 0.5 and \
            (sp != "neoformans" or r.contig.iloc[0] == "NC_026751.1")
        if not ok:
            b.add_patch(Rectangle((j - 0.45, i - 0.45), 0.9, 0.9, fc="white", ec="#D0D0D0", lw=0.5))
            b.text(j, i, "–", ha="center", va="center", fontsize=7, color=LOST)
            continue
        v = 100 * r.identity.iloc[0]
        b.add_patch(Rectangle((j - 0.45, i - 0.45), 0.9, 0.9, fc=cmap(norm(v)), ec="white", lw=0.5))
        star = "*" if r.broken.iloc[0] else ""
        b.text(j, i, f"{v:.0f}{star}", ha="center", va="center", fontsize=5.5,
               color="white" if v > 92 else SC_AXIS)
b.set_xlim(-0.5, len(cols) - 0.5)
b.set_ylim(len(BLOCK) - 0.5, -0.5)
b.set_xticks(range(len(cols)), [n for n, _ in cols], fontsize=5.3, rotation=0)
b.xaxis.set_ticks_position("top")
for lab, (n, sp) in zip(b.get_xticklabels(), cols):
    if n.startswith("C."):
        lab.set_style("italic")
DESC = {"CGB_F4640C": "4640 (AA transporter)", "CGB_F4660C": "4660 (dihydropyrimidinase)",
        "CGB_F4610W": "HOC1 (mannosyltransf.)", "CGB_F4670W": "4670 (2OG oxygenase)", "CGB_F4690W": "RGD1 (GAP)"}
b.set_yticks(range(len(BLOCK)), [DESC.get(g, SHORT[g]) for g in BLOCK], fontsize=6)
for lab, g in zip(b.get_yticklabels(), BLOCK):
    if g in LOST_IN_NEO:
        lab.set_color(LOST)
b.tick_params(length=0)
for s in b.spines.values():
    s.set_visible(False)
panel(b, "B", dx=-0.2, dy=1.13)

# C: block tree
tree = Phylo.read("vnii/tree/block.contree", "newick")
tree.root_with_outgroup("deneoformans_ref")
out = max(tree.root.clades, key=lambda k: k.branch_length or 0)  # the long root branch, drawn short
true_len, out.branch_length = out.branch_length, 0.004
tips = tree.get_terminals()
y = {t: i for i, t in enumerate(tips)}
depth = tree.depths()


def ypos(cl):
    if cl not in y:
        y[cl] = sum(ypos(k) for k in cl.clades) / len(cl.clades)
    return y[cl]


ypos(tree.root)
for cl in tree.find_clades(order="preorder"):
    for k in cl.clades:
        c.plot([depth[cl], depth[cl]], [y[cl], y[k]], color=SC_AXIS, lw=0.6)
        c.plot([depth[cl], depth[k]], [y[k], y[k]], color=SC_AXIS, lw=0.6)
        if not k.is_terminal() and k.confidence is not None and k.confidence >= 95 and len(k.get_terminals()) >= 5:
            c.text(depth[k] - 0.0004, y[k] - 0.35, f"{k.confidence:.0f}", fontsize=5, ha="right", va="bottom",
                   color=SC_ANNOT)


# VGI tips are named from the sample table by run, not from the tip name baked into the tree, so a
# corrected strain name (census/fix_strain_names.py) reaches the figure without rebuilding the tree
VGI_NAME = {r.run: (r.country, r.strain) for r in pd.read_csv("vnii/donor/vgi_samples.tsv", sep="\t").itertuples()}
SHORT_COUNTRY = {"Papua New Guinea": "PNG", "South Africa": "S. Africa"}


def tip_label(n):
    if n.startswith("VNII_"):
        return "VNII " + n[5:], LOST
    if n.startswith("VGI_"):
        run = n.rsplit("_", 1)[1]
        if run in VGI_NAME:
            country, strain = VGI_NAME[run]
            return f"VGI {SHORT_COUNTRY.get(country, country)} {strain}".replace("_", " "), GATTII
        return "VGI " + n[4:].replace("_ref", "").replace("_", " "), GATTII  # the WM276 reference
    return "C. " + n.replace("_ref", ""), SC_AXIS


xmax = max(depth[t] for t in tips)
for t in tips:
    lab, colour = tip_label(t.name)
    c.text(depth[t] + xmax * 0.02, y[t], lab, va="center", fontsize=5.5, color=colour,
           style="italic" if lab.startswith("C. ") else "normal")
c.set_xlim(0, xmax * 1.75)
c.set_ylim(len(tips) + 1.6, -1)
c.set_yticks([])
c.spines["left"].set_visible(False)
ob = ypos(out)
c.text(depth[out] * 0.5, ob, "//", ha="center", va="center", fontsize=7, color=SC_AXIS,
       bbox=dict(fc="white", ec="none", pad=0.2))
c.text(depth[out] * 0.5, ob + 0.9, f"{true_len:.3f}", ha="center", va="top", fontsize=5, color=SC_ANNOT)
c.plot([xmax * 0.95, xmax * 0.95 + 0.002], [len(tips) + 0.3] * 2, color=SC_AXIS, lw=0.8)
c.text(xmax * 0.95 + 0.001, len(tips) + 0.5, "0.002 subst./site", ha="center", va="top", fontsize=5.5)
c.set_xticks([])
c.spines["bottom"].set_visible(False)
panel(c, "C", dx=-0.02, dy=1.0)


# D: divergence along the locus (added 30 Oct 2026; data in figures/msystems/data/)
REVD = "figures/msystems/data/"
wins = pd.read_csv(REVD + "identity_windows.tsv", sep="\t")
wins["mid"] = (wins.start + 500) / 1000
gdiv = pd.read_csv(REVD + "gene_divergence.tsv", sep="\t")
gpos = {}
for line in open(REVD + "genes_vs_vnii_loc.paf"):
    f = line.split("\t")
    gpos.setdefault(f[0], ((int(f[7]) + int(f[8])) / 2000, int(f[7]) / 1000, int(f[8]) / 1000))
d.axvspan(79.9, 90.8, color=KDZ, alpha=0.30, lw=0)
d.axvspan(90.8, 114.9, color=GATTII, alpha=0.10, lw=0)
for col, colour, lab in (("h99", NEO, "H99 (VNI)"), ("wm", GATTII, "WM276 (VGI)")):
    w = wins[wins["cov_" + col] >= 0.5]
    d.plot(w.mid, 100 * (1 - w["id_" + col]), "o", ms=2.2, color=colour, lw=0, label="VNII vs " + lab)
d.axhline(1.58, color=SC_ANNOT, ls=":", lw=0.7)
d.text(60, 5.4, "VNII\u2013H99 genome-wide 1.6%", fontsize=5.5, ha="center", va="bottom", color=SC_ANNOT)
g = gdiv.dropna(subset=["VNII-H99 %"])
d.scatter([gpos[t][0] for t in g.wm_tag], g["VNII-H99 %"], marker="D", s=20, facecolor="white",
          edgecolor=NEO, lw=0.9, zorder=5, label="VNII vs H99, whole gene")
for tag in BLOCK + ["CGB_F4600C", "CGB_F4700C"]:
    if tag in gpos:
        mid, s0, s1 = gpos[tag]
        inb = tag in BLOCK
        d.add_patch(Rectangle((s0, -3.1), s1 - s0, 1.5, color=GATTII if inb else FLANK, lw=0))
for tag, nm, ha, off in (("CGB_F4600C", "F4600C", "center", 0), ("CGB_F4610W", "HOC1", "center", 0),
                         ("CGB_F4690W", "RGD1", "right", -0.8), ("CGB_F4700C", "F4700C", "left", 1.2)):
    d.text(gpos[tag][0] + off, -5.0, nm, ha=ha, fontsize=5.5, color=SC_ANNOT)
d.text(85.3, 22.0, "KDZ\ntransposon", ha="center", fontsize=5.5, color=SC_ANNOT)
d.text(102.9, 23.2, "$\\it{gattii}$-type block", ha="center", fontsize=5.5, color=GATTII)
d.set_xlim(12, 183)
d.set_ylim(-6.0, 26)
d.set_xlabel("Position along the VNII locus (kb; long-read genome GCA_022832995.1, chr. 5)")
d.set_ylabel("Substitutions\nper site (%)")
d.legend(loc="upper left", fontsize=5.5, handletextpad=0.3, borderpad=0.2, labelspacing=0.25)
panel(d, "D", dx=0.0, dy=1.02)

for ext, dpi in (("pdf", 600), ("png", 300), ("tiff", 300)):
    fig.savefig(f"figures/msystems/out/Figure5.{ext}", dpi=dpi)
print("written: figures/msystems/out/Figure5")
