"""3'-anchored (MAMA-style) primer design from cladeprint markers.

Marker regions after near-miss filtering are short (40-44 bp on average), too short to hold a
conventional primer pair. Discrimination in allele-specific PCR comes from the 3' end, so instead
of putting both primers inside marker sequence this designs one primer whose 3' end sits on
group-specific bases and lets its partner sit in conserved flanking sequence, as the published
MAMA assays for the C. gattii complex do (Kelley et al. 2014).

  windows <bundle> <samples.tsv> <genome_dir> <out>   marker regions + flanks from a member genome
  design  <out>                                       Primer3 pairs anchored in the marker region

Step `windows` locates each marker unitig in one representative member genome by its first 31
bases (via examples/primer_sites, which also reports how many equally good sites there are, so
multi-copy regions are dropped) and cuts a window of FLANK bases either side.
Step `design` runs Primer3 on each window with the product forced to end inside the marker region,
so the anchored primer's 3' end lands on group-specific sequence.
"""

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd
import primer3

FLANK = 220
MIN_MARKER = 34          # need a few positions to anchor on
TOP_UNITIGS = 150        # longest marker unitigs per group
N_RETURN = 5
PRODUCT = (90, 400)
PRIMER_SITES = str(Path.home() / "cladeprint/target/release/examples/primer_sites")
LEVELS = ["species", "lineage", "sublineage", "subclade"]


def read_fasta(path):
    name, seq = None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                yield name, "".join(seq)
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
    if name:
        yield name, "".join(seq)


def group_paths(s):
    out = {}
    for lv in LEVELS:
        if lv not in s.columns:
            continue
        cols = LEVELS[: LEVELS.index(lv) + 1]
        p = s[cols].apply(lambda r: None if any(pd.isna(x) or x == "" for x in r) else "/".join(r), axis=1)
        for g in p.dropna().unique():
            out[g] = list(p.index[p == g])
    return out


def windows(args):
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    s = pd.read_csv(args.samples, sep="\t", dtype=str).set_index("sample")
    qc = pd.read_csv(Path(args.samples).parent / "qc.tsv", sep="\t", dtype=str).set_index("sample")
    s = s[qc.reindex(s.index).train.values == "true"]
    members = group_paths(s)
    peak = pd.to_numeric(qc.peak, errors="coerce")
    rows = []
    for fa in sorted(Path(args.bundle, "markers").glob("*/*.fa")):
        group = fa.stem.replace("__", "/")
        if group not in members:
            print(f"[skip] {group}: no members in the sheet", file=sys.stderr)
            continue
        rep = max(members[group], key=lambda x: (peak.get(x, 0), x))
        unis = [(n, q) for n, q in read_fasta(fa) if len(q) >= MIN_MARKER and set(q) <= set("ACGT")]
        unis.sort(key=lambda x: -len(x[1]))
        unis = unis[: args.top]
        if not unis:
            print(f"[skip] {group}: no marker unitig >= {MIN_MARKER} bp", file=sys.stderr)
            continue
        rows.append((group, rep, unis))
    # one primer_sites call per group against its representative genome
    recs = []
    for group, rep, unis in rows:
        tag = group.replace("/", "_")
        od = out / "loc" / tag
        od.mkdir(parents=True, exist_ok=True)
        with open(od / "oligos.tsv", "w") as f:
            f.write("oligo\tsequence\n")
            for n, q in unis:
                f.write(f"{n}\t{q[:31]}\n")
        with open(od / "genome.tsv", "w") as f:
            f.write(f"{rep}\t{Path(args.genome_dir) / (rep + '.fa')}\n")
        subprocess.run([PRIMER_SITES, str(od / "oligos.tsv"), str(od / "genome.tsv"), str(od / "sites.tsv")],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        sites = pd.read_csv(od / "sites.tsv", sep="\t", dtype={"mismatches": str}).set_index("oligo")
        gseq = dict(read_fasta(Path(args.genome_dir) / (rep + ".fa")))
        useq = dict(unis)
        n_ok = 0
        with open(out / f"{tag}.windows.fa", "w") as f:
            for n, _ in unis:
                if n not in sites.index:
                    continue
                r = sites.loc[n]
                if str(r.mismatches) != "0" or int(r.n_sites) != 1:
                    continue  # absent, imperfect or multi-copy in the representative genome
                contig, pos = r.contig, int(r.pos)
                c = gseq.get(contig)
                if c is None:
                    continue
                mlen = len(useq[n])
                if r.strand == "-":
                    continue  # keep it simple: marker found on the + strand of this contig
                a, b = max(0, pos - FLANK), min(len(c), pos + mlen + FLANK)
                win = c[a:b]
                mstart, mend = pos - a, pos - a + mlen
                if mend > len(win) or mstart < 0:
                    continue
                f.write(f">{tag}|{n}|rep={rep}|marker={mstart}-{mend}|winlen={len(win)}\n{win}\n")
                n_ok += 1
        recs.append(dict(group=group, level=fa.parent.name, rep=rep, unitigs=len(unis), windows=n_ok))
        print(f"{group}: rep {rep}, {n_ok}/{len(unis)} windows")
    pd.DataFrame(recs).to_csv(out / "windows_summary.tsv", sep="\t", index=False)


def design(args):
    """Force the right primer's 3' end onto group-specific bases.

    Primer3 reports a right primer by its 5' end, so the 3' end is pinned with
    SEQUENCE_FORCE_RIGHT_END instead of inferred. Up to ANCHORS positions spread across each
    marker region are tried; the left primer is free to sit in conserved flanking sequence.
    """
    out = Path(args.out)
    rows = []
    for fa in sorted(out.glob("*.windows.fa")):
        for head, win in read_fasta(fa):
            parts = dict(p.split("=", 1) for p in head.split("|") if "=" in p)
            tag, uni = head.split("|")[0], head.split("|")[1]
            ms, me = (int(x) for x in parts["marker"].split("-"))
            span = max(1, me - ms)
            step = max(1, span // args.anchors)
            for anchor in range(ms, me, step):
                try:
                    res = primer3.bindings.design_primers(
                        {"SEQUENCE_ID": head, "SEQUENCE_TEMPLATE": win,
                         "SEQUENCE_FORCE_RIGHT_END": anchor},
                        {"PRIMER_TASK": "generic", "PRIMER_PICK_LEFT_PRIMER": 1,
                         "PRIMER_PICK_RIGHT_PRIMER": 1, "PRIMER_NUM_RETURN": 1,
                         "PRIMER_PRODUCT_SIZE_RANGE": [list(PRODUCT)]},
                    )
                except Exception as e:
                    print(f"[warn] {head} @{anchor}: {e}", file=sys.stderr)
                    continue
                if not res.get("PRIMER_PAIR_NUM_RETURNED", 0):
                    continue
                rstart, rlen = res["PRIMER_RIGHT_0"]
                rows.append(dict(
                    group=tag, window=head, unitig=uni,
                    pair=f"{tag}|{uni}|{anchor}",
                    left=res["PRIMER_LEFT_0_SEQUENCE"], right=res["PRIMER_RIGHT_0_SEQUENCE"],
                    left_tm=round(res["PRIMER_LEFT_0_TM"], 1), right_tm=round(res["PRIMER_RIGHT_0_TM"], 1),
                    product_bp=res["PRIMER_PAIR_0_PRODUCT_SIZE"],
                    anchor_pos=anchor, right_3p=rstart - rlen + 1,
                    marker_start=ms, marker_end=me,
                    anchored=ms <= anchor < me,
                    penalty=round(res["PRIMER_PAIR_0_PENALTY"], 3),
                ))
    d = pd.DataFrame(rows)
    if len(d):
        d = d.sort_values("penalty").drop_duplicates(["group", "right"])
    d.to_csv(out / "pairs_3p.tsv", sep="\t", index=False)
    print(f"{len(d)} pairs, {int(d.anchored.sum()) if len(d) else 0} anchored in the marker region")
    if len(d):
        print(d[d.anchored].groupby("group").size().to_string())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    w = sub.add_parser("windows")
    w.add_argument("bundle"); w.add_argument("samples"); w.add_argument("genome_dir"); w.add_argument("out")
    w.add_argument("--top", type=int, default=TOP_UNITIGS); w.set_defaults(f=windows)
    d = sub.add_parser("design"); d.add_argument("out")
    d.add_argument("--anchors", type=int, default=6); d.set_defaults(f=design)
    a = ap.parse_args()
    a.f(a)
