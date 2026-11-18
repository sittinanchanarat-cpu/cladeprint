# helper for Figure 6A and Figure S6 of the paper.
"""Assay figure, panel A: published assays x genome groups, in-silico PCR on NCBI assemblies.

Cell = share of the group's assemblies giving the expected product (lenient 3' rule; probes
must sit in the product). A dot marks cells where the strict rule gives at least 0.1 less.
Boxed = intended target. Multiplex typing assays print their commonest band sizes (bp).
Columns with fewer than 5 assemblies are greyed in the header.
MAIN_GROUPS picks the main-text rows; --all renders every group for the supplement.
ARTIFACT marks cells our 3' rule scores 0 although the assay is known to work (URA5 in VGII:
a mismatch at the penultimate base of the forward primer, which PCR tolerates).
Input: census/asm/r5_pcr_products.tsv, classify_subclade.tsv, assemblies.tsv.
Usage: python figures/msystems/assay_panel.py [--all] -> figures/msystems/out/assay_panel[_full].{pdf,png,tiff}
"""
import sys
from collections import Counter

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SC_AXIS, SC_ANNOT = "#2C2C2C", "#4A4A4A"
SC_WARM = mpl.colors.LinearSegmentedColormap.from_list("sc_warm", ["#F5ECD7", "#DFB166", "#C96A3A", "#9B2335"])
mpl.rcParams.update({
    "font.family": ["Arial", "Liberation Sans", "DejaVu Sans"], "font.size": 7,
    "axes.linewidth": 0.75, "axes.edgecolor": SC_AXIS, "xtick.color": SC_AXIS, "ytick.color": SC_AXIS,
    "xtick.major.size": 0, "ytick.major.size": 0, "figure.facecolor": "white", "savefig.facecolor": "white",
    "pdf.fonttype": 42,
    # mathtext ($\it{}$ species names) defaults to DejaVu: force the same family as the rest
    "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
})

ASM = "census/asm/"
P = "Lengeler 2001 serotype/mating-type allele panel: "
Y = "Yan 2002 STE20/STE12 serotype + mating-type PCR: "
# (reaction, row label, intended columns, pattern assay?)
SPECIES = ["C.neoformans"]
ROWS = [
    ("Species", [
        ("Aoki 1999 serotype A RAPD-derived PCR (CNa-70)", "Aoki CNa-70 (serotype A)", "Cn", False),
        ("Aoki 1999 serotype A RAPD-derived PCR (CNa-29)", "Aoki CNa-29 (serotype A)", "Cn", False),
        ("Aoki 1999 serotype B PCR (CNb-49)", "Aoki CNb-49 (serotype B)", "gattii", False),
        ("Mbangiwa 2024 QSP1 qPCR (C. neoformans)", "QSP1 qPCR, @C. neoformans@", "Cn", False),
        ("Mbangiwa 2024 QSP1 qPCR (C. deneoformans)", "QSP1 qPCR, @C. deneoformans@", "Cd", False),
        ("Mbangiwa 2024 QSP1 qPCR (C. gattii)", "QSP1 qPCR, @C. gattii@", "gattii", False),
    ]),
    ("Serotyping\n(band size, bp)", [
        ("Feng 2013 STR1 singleplex species PCR", "Feng STR1", "all", True),
        ("Enache-Angoulvant 2007 CAP59 PCR-REA / fragment-length serotyping", "Enache CAP59", "all", True),
        ("Halliday 1999 MFalpha/660 mating-type + variety coamplification", "Halliday MFα / 660", "all", True),
    ]),
    ("Serotype alleles", [
        (P + "A053+A052", "CLA4-A (Lengeler)", "Cn", False),
        (P + "JOHE3065+JOHE3066", "CLA4-D (Lengeler)", "Cd", False),
        (P + "A054+A056", "CNA1-A (Lengeler)", "Cn", False),
        (P + "JOHE2926+JOHE3238", "CNA1-D (Lengeler)", "Cd", False),
        (P + "A057+A059", "GPA1-A (Lengeler)", "Cn", False),
        (P + "JOHE2596+JOHE3240", "GPA1-D (Lengeler)", "Cd", False),
    ]),
    ("Serotype × MAT", [
        (Y + "STE20Aalpha Forward+STE20Aalpha Reverse", "STE20α-A (Yan)", "Cn", False),
        (Y + "STE20Aa Forward+STE20Aa Reverse", "STE20a-A (Yan)", "Cn", False),
        (Y + "STE20Dalpha Forward+STE20Dalpha Reverse", "STE20α-D (Yan)", "Cd", False),
        (Y + "STE20Da Forward+STE20Da Reverse", "STE20a-D (Yan)", "Cd", False),
        (P + "JOHE1909+JOHE1910", "MFα2-D (Lengeler)", "Cd", False),
    ]),
    ("Molecular typing", [
        ("URA5-RFLP (Meyer et al. 2003)", "URA5-RFLP amplicon", "all", False),
        ("PLB1-RFLP (Latouche et al. 2003)", "PLB1-RFLP amplicon", "all", False),
    ]),
]
MAIN_GROUPS = ("Species", "Serotype alleles")
ARTIFACT = {("URA5-RFLP (Meyer et al. 2003)", "Cdg")}
COLS = [("VNI", "VNI"), ("VNII", "VNII"), ("VNB", "VNB"), ("Cd", "@C. deneo-@\n@formans@"),
        ("Cg", "@C. gattii@"), ("Cdg", "@C. deutero-@\n@gattii@"), ("Cb", "@C. bacilli-@\n@sporus@"),
        ("Ct", "@C. tetra-@\n@gattii@"), ("Cdc", "@C. deca-@\n@gattii@")]
SP = {"C.deneoformans": "Cd", "C.gattii": "Cg", "C.deuterogattii": "Cdg", "C.bacillisporus": "Cb",
      "C.tetragattii": "Ct", "C.decagattii": "Cdc"}
TARGET = {"Cn": {"VNI", "VNII", "VNB"}, "Cd": {"Cd"}, "gattii": {"Cg", "Cdg", "Cb", "Ct", "Cdc"},
          "all": {c for c, _ in COLS if c != "AD"}}


def ital(s):
    """@species name@ -> mathtext italics, one segment per line (mathtext cannot span newlines)."""
    out = []
    for i, part in enumerate(s.split("@")):
        if i % 2 == 0 or not part:
            out.append(part)
        else:
            tail = "-" if part.endswith("-") else ""
            out.append(r"$\it{" + part[: len(part) - len(tail)].replace(" ", r"\,") + "}$" + tail)
    return "".join(out)


def column_of(g, cl, meta):
    if " x " in meta.at[g, "organism"]:
        return "AD"
    fl = cl.at[g, "flags"]
    if "EXCESS_KMERS" in fl and "MIXED" in fl:
        return "AD"   # hybrid by k-mer excess + mixed parents, whatever NCBI calls it
    sp = cl.at[g, "species"]
    if sp == "C.neoformans":
        lin = cl.at[g, "lineage"].split("/")[-1]
        return lin if lin in ("VNI", "VNII", "VNB") else None
    return SP.get(sp)


def main(full=False, ax=None, stem=None, label="A"):
    """Standalone figure, or with ax: draw into that axis (Fig. 6A) and return the row extent.
    stem/label: output path stem and panel letter of the standalone figure (Fig. S6: no letter)."""
    groups = ROWS if full else [(g, rs) for g, rs in ROWS if g in MAIN_GROUPS]
    prod = pd.read_csv(ASM + "r5_pcr_products.tsv", sep="\t")
    prod = prod[prod.expected & prod.probe_ok]
    cl = pd.read_csv(ASM + "classify_subclade.tsv", sep="\t", dtype=str, keep_default_na=False).set_index("sample")
    meta = pd.read_csv(ASM + "assemblies.tsv", sep="\t").set_index("accession")
    col = {g: column_of(g, cl, meta) for g in cl.index}
    members = {c: [g for g, x in col.items() if x == c] for c, _ in COLS}
    n = {c: len(v) for c, v in members.items()}

    rows = [(grp, *r) for grp, rs in groups for r in rs]
    val = np.full((len(rows), len(COLS)), np.nan)
    strict = np.full_like(val, np.nan)
    sizes = {}
    for i, (_, rx, _, _, pattern) in enumerate(rows):
        for rule, arr in (("lenient", val), ("strict", strict)):
            have = set(prod[(prod.reaction == rx) & (prod.rule == rule)].genome)
            for j, (c, _) in enumerate(COLS):
                if n[c]:
                    arr[i, j] = len(have & set(members[c])) / n[c]
        if pattern:
            p = prod[(prod.reaction == rx) & (prod.rule == "lenient")]
            for j, (c, _) in enumerate(COLS):
                pats = Counter("/".join(map(str, sorted({int(round(s, -1)) for s in x})))
                               for _, x in p[p.genome.isin(members[c])].groupby("genome")["size"])
                sizes[i, j] = pats.most_common(1)[0][0] if pats else ""

    # layout: one gap row between groups
    ypos, y, ticks = [], 0, []
    for k, (grp, rs) in enumerate(groups):
        if k:
            y += 0.6
        start = y
        for _ in rs:
            ypos.append(y)
            y += 1
        ticks.append((grp, start, y - 1))
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(7.0, 0.19 * y + 1.3))
    else:
        fig = ax.figure
    for i, yy in enumerate(ypos):
        for j in range(len(COLS)):
            v = val[i, j]
            face = "white" if np.isnan(v) else SC_WARM(v)
            ax.add_patch(plt.Rectangle((j, yy), 1, 1, facecolor=face, edgecolor="white", lw=1.2))
            if np.isnan(v):
                continue
            if (rows[i][1], COLS[j][0]) in ARTIFACT:
                ax.add_patch(plt.Rectangle((j, yy), 1, 1, facecolor="none", edgecolor="#BFBFBF",
                                           hatch="////", lw=0))
                ax.text(j + 0.5, yy + 0.5, "†", ha="center", va="center", fontsize=7, color=SC_ANNOT)
                continue
            txt = sizes.get((i, j)) if rows[i][4] else f"{v:.2f}".lstrip("0") if v < 1 else "1"
            if rows[i][4] and not txt:
                txt = "–"
            ax.text(j + 0.5, yy + 0.5, txt, ha="center", va="center", fontsize=5.5 if rows[i][4] else 6,
                    color="white" if v > 0.62 else SC_AXIS)
            if strict[i, j] < v - 0.1:
                ax.plot(j + 0.88, yy + 0.18, "o", ms=1.8, color=SC_ANNOT, mec="white", mew=0.3)
        # intended target outline (contiguous runs)
        tg = TARGET[rows[i][3]]
        js = [j for j, (c, _) in enumerate(COLS) if c in tg]
        for a, b in runs(js):
            ax.add_patch(plt.Rectangle((a + 0.04, ypos[i] + 0.06), b - a + 0.92, 0.88, fill=False,
                                       edgecolor=SC_AXIS, lw=0.9))
    ax.set_xlim(0, len(COLS))
    ax.set_ylim(y, -0.1)
    ax.set_yticks([yy + 0.5 for yy in ypos])
    ax.set_yticklabels([ital(r[2]) for r in rows])
    ax.set_xticks([j + 0.5 for j in range(len(COLS))])
    ax.set_xticklabels([ital(lab) + f"\n(n={n[c]})" for c, lab in COLS], fontsize=6.2)
    ax.xaxis.tick_top()
    for t, (c, _) in zip(ax.get_xticklabels(), COLS):
        if n[c] < 5:
            t.set_color("#9A9A9A")
    ax.spines[:].set_visible(False)
    # brace for C. neoformans lineages
    ax.annotate("", xy=(0.05, -2.55), xytext=(2.95, -2.55), xycoords="data", annotation_clip=False,
                arrowprops=dict(arrowstyle="-", color=SC_AXIS, lw=0.75))
    ax.text(1.5, -2.75, "C. neoformans", ha="center", va="bottom", fontstyle="italic", fontsize=6.5)
    # group labels at the right
    for grp, a, b in ticks:
        ax.plot([len(COLS) + 0.25] * 2, [a + 0.15, b + 0.85], color=SC_AXIS, lw=0.75, clip_on=False)
        ax.text(len(COLS) + 0.4, (a + b + 1) / 2, grp, ha="left", va="center", fontsize=6.5)
    cax = fig.add_axes([0.80, 0.06, 0.12, 0.012]) if own else ax.inset_axes([0.72, -0.06, 0.22, 0.012])
    cb = fig.colorbar(mpl.cm.ScalarMappable(cmap=SC_WARM, norm=mpl.colors.Normalize(0, 1)), cax=cax,
                      orientation="horizontal", ticks=[0, 0.5, 1])
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=6, length=2)
    cb.set_label("share with expected product", fontsize=6.5)
    if not own:
        return y
    if label:
        fig.text(0.0, 0.985, label, fontsize=10, fontweight="bold", va="top")
    stem = stem or ("figures/msystems/out/assay_panel_full" if full else "figures/msystems/out/assay_panel")
    fig.savefig(stem + ".pdf", dpi=600, bbox_inches="tight")
    fig.savefig(stem + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(stem + ".tiff", dpi=300, bbox_inches="tight")
    print({c: n[c] for c, _ in COLS})


def runs(js):
    out, s = [], None
    for k, j in enumerate(js):
        if s is None:
            s = j
        if k == len(js) - 1 or js[k + 1] != j + 1:
            out.append((s, j))
            s = None
    return out


if __name__ == "__main__":
    main(full="--all" in sys.argv)
