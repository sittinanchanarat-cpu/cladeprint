"""Primer pairs from cladeprint markers: Primer3 for the thermodynamics, a census-wide specificity
screen for the part Primer3 cannot do.

Both primers of a pair are taken from inside one marker unitig, so every base is group-specific
sequence and the product size is fixed by construction. Specificity is then judged against every
census genome that is not a member of the group: a pair is kept only if each primer's *best* site
in every non-member genome is at least MIN_MM mismatches away, with at least MIN_MM_3P of them in
the 3' TAIL bases. No non-member can give a product if neither primer binds anywhere in it.

  design <bundle_dir> <out_dir> [--per-group N]
      -> out_dir/pairs.tsv (candidate pairs) and out_dir/oligos.tsv (for examples/primer_sites)
  examples/primer_sites out_dir/oligos.tsv genomes.tsv out_dir/sites.tsv
  screen <out_dir> <samples.tsv> [--level-cols species,lineage,sublineage,subclade]
      -> out_dir/pairs_screened.tsv (per pair: worst non-member distance, members detected, verdict)
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd
import primer3

MIN_UNITIG = 120
PRODUCT = (90, 300)
PER_GROUP = 40
N_RETURN = 3
MIN_MM = 3
MIN_MM_3P = 1
TAIL = 5  # must match examples/primer_sites


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


def design(args):
    bundle, out = Path(args.bundle), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for fa in sorted((bundle / "markers").glob("*/*.fa")):
        group = fa.stem.replace("__", "/")
        level = fa.parent.name
        cands = [(n, s) for n, s in read_fasta(fa) if len(s) >= MIN_UNITIG and re.fullmatch(r"[ACGT]+", s)]
        cands.sort(key=lambda x: -len(x[1]))
        n_pairs = 0
        for name, seq in cands:
            if n_pairs >= args.per_group:
                break
            try:
                res = primer3.bindings.design_primers(
                    {"SEQUENCE_ID": name, "SEQUENCE_TEMPLATE": seq},
                    {
                        "PRIMER_TASK": "generic",
                        "PRIMER_PICK_LEFT_PRIMER": 1,
                        "PRIMER_PICK_RIGHT_PRIMER": 1,
                        "PRIMER_NUM_RETURN": N_RETURN,
                        "PRIMER_PRODUCT_SIZE_RANGE": [list(PRODUCT)],
                    },
                )
            except Exception as e:  # a template Primer3 cannot use
                print(f"[warn] {group} {name}: {e}", file=sys.stderr)
                continue
            for i in range(res.get("PRIMER_PAIR_NUM_RETURNED", 0)):
                if n_pairs >= args.per_group:
                    break
                rows.append(
                    dict(
                        group=group,
                        level=level,
                        unitig=name,
                        unitig_len=len(seq),
                        pair=f"{group}|{name}|{i}".replace("/", "_"),
                        left=res[f"PRIMER_LEFT_{i}_SEQUENCE"],
                        right=res[f"PRIMER_RIGHT_{i}_SEQUENCE"],
                        left_tm=round(res[f"PRIMER_LEFT_{i}_TM"], 1),
                        right_tm=round(res[f"PRIMER_RIGHT_{i}_TM"], 1),
                        left_gc=round(res[f"PRIMER_LEFT_{i}_GC_PERCENT"], 1),
                        right_gc=round(res[f"PRIMER_RIGHT_{i}_GC_PERCENT"], 1),
                        product_bp=res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
                        penalty=round(res[f"PRIMER_PAIR_{i}_PENALTY"], 3),
                    )
                )
                n_pairs += 1
    p = pd.DataFrame(rows)
    p.to_csv(out / "pairs.tsv", sep="\t", index=False)
    olig = pd.concat(
        [
            p[["pair", "left"]].rename(columns={"pair": "oligo", "left": "sequence"}).assign(oligo=lambda d: d.oligo + "|L"),
            p[["pair", "right"]].rename(columns={"pair": "oligo", "right": "sequence"}).assign(oligo=lambda d: d.oligo + "|R"),
        ]
    ).drop_duplicates("sequence")
    olig.to_csv(out / "oligos.tsv", sep="\t", index=False)
    print(f"{len(p)} pairs over {p.group.nunique()} groups; {len(olig)} unique primers -> {out}")
    print(p.groupby(["level", "group"]).size().to_string())


def screen(args):
    out = Path(args.out)
    pairs = pd.read_csv(out / "pairs.tsv", sep="\t")
    sites = pd.read_csv(out / "sites.tsv", sep="\t", dtype={"mismatches": str})
    s = pd.read_csv(args.samples, sep="\t", dtype=str)
    levels = args.level_cols.split(",")
    s = s.set_index("sample")
    seq_of = {}
    for _, r in pairs.iterrows():
        seq_of[(r.pair, "L")] = r.left
        seq_of[(r.pair, "R")] = r.right
    # sites are keyed by sequence (duplicates collapsed): map back by sequence
    olig = pd.read_csv(out / "oligos.tsv", sep="\t")
    by_seq = dict(zip(olig.oligo, olig.sequence))
    sites["sequence"] = sites.oligo.map(by_seq)
    sites["mm"] = pd.to_numeric(sites.mismatches, errors="coerce").fillna(99).astype(int)
    sites["tail"] = pd.to_numeric(sites.mm_3prime, errors="coerce").fillna(99).astype(int)
    key = sites.set_index(["sequence", "sample"])[["mm", "tail"]]

    def path(row, level):
        cols = levels[: levels.index(level) + 1]
        v = [row.get(c) for c in cols]
        return None if any(pd.isna(x) or x in (None, "") for x in v) else "/".join(v)

    member = {}
    for lv in levels:
        p = s.apply(lambda r: path(r, lv), axis=1)
        for g in p.dropna().unique():
            member[g] = set(p.index[p == g])
    rows = []
    for _, r in pairs.iterrows():
        mem = member.get(r.group, set())
        non = [x for x in s.index if x not in mem]
        rec = {}
        for side, seq in (("L", r.left), ("R", r.right)):
            sub = key.loc[seq] if seq in key.index.get_level_values(0) else None
            if sub is None:
                rec[side] = (99, 99, 0.0)
                continue
            non_mm = sub.reindex(non).mm.fillna(99)
            non_tail = sub.reindex(non).tail.fillna(99)
            mem_mm = sub.reindex(list(mem)).mm.fillna(99)
            rec[side] = (int(non_mm.min()), int(non_tail[non_mm == non_mm.min()].min()), float((mem_mm == 0).mean()))
        worst = min(rec["L"][0], rec["R"][0])
        tail_at_worst = min(rec["L"][1], rec["R"][1])
        detect = min(rec["L"][2], rec["R"][2])
        rows.append(
            dict(
                pair=r.pair, group=r.group, level=r.level, product_bp=r.product_bp,
                left=r.left, right=r.right, left_tm=r.left_tm, right_tm=r.right_tm,
                n_members=len(mem), member_detect=round(detect, 4),
                nonmember_min_mm=worst, nonmember_mm_3prime=tail_at_worst,
                pass_specificity=bool(worst >= MIN_MM and tail_at_worst >= MIN_MM_3P),
            )
        )
    d = pd.DataFrame(rows).sort_values(["group", "pass_specificity", "member_detect"], ascending=[True, False, False])
    d.to_csv(out / "pairs_screened.tsv", sep="\t", index=False)
    print(d.groupby("group").agg(pairs=("pair", "size"), passing=("pass_specificity", "sum"),
                                 best_detect=("member_detect", "max")).to_string())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("design"); d.add_argument("bundle"); d.add_argument("out")
    d.add_argument("--per-group", type=int, default=PER_GROUP); d.set_defaults(f=design)
    g = sub.add_parser("screen"); g.add_argument("out"); g.add_argument("samples")
    g.add_argument("--level-cols", default="species,lineage,sublineage,subclade"); g.set_defaults(f=screen)
    a = ap.parse_args()
    a.f(a)
