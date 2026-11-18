"""Table 2: leave-one-clade-out (R2) calls under the current classify rules.

The R2 runs (18 Oct) were classified before the rule changes of 19 Oct (high confidence from the
LOO 5th percentile; a low call stops the descent unless LOW_COVERAGE) and 22 Oct (at LOW_COVERAGE
a call is no more confident than its parent). The per-group scores they stored do not depend on
the rules, so the calls are re-made here from those scores with a re-implementation of
classify.rs::call. With --check the old rules are applied instead and every stored call and
confidence must be reproduced, which tests the re-implementation.

Runs counted: the held-out runs that pass QC in the census run (as in r2_summary.py), literature
labels and second-round labels together. Right answer at the held-out level: no call for a whole
held-out lineage or sublineage (NO_CALL, or NOT_DESCENDED below a low call); the still-trained
parent group for a VNIa sub-clade or a deep clade of a species ("other1", "other2").

Usage: python census/r2_table.py census/out [--check]
Writes census/out/r2/table2_long.tsv (one row per holdout x level) and prints Table 2.
"""

import csv
import sys
import tomllib
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median

CALL_MIN = MIXED_MIN = 0.3
LEVELS = ["species", "lineage", "sublineage"]
RANK = {"low": 0, "medium": 1, "high": 2}


def related(a, b):
    return a == b or a.startswith(b + "/") or b.startswith(a + "/")


def p05(v):
    v = sorted(v)
    return v[round((len(v) - 1) * 0.05)]


def call(groups, levels, sc, low_cov, old=False):
    """classify.rs::call on stored scores. groups: id -> {level, loo, p05}; sc: id -> score."""
    out, parent, parent_conf, stopped = [], None, None, False
    for lv, _ in enumerate(levels):
        at = [(g, s) for g, s in sc.items() if groups[g]["level"] == lv and s is not None]
        cand = [(g, s) for g, s in at if (lv == 0 if parent is None else g.startswith(parent + "/"))]
        best = max(cand, key=lambda x: x[1]) if cand else None
        if stopped:
            out.append((None, "none", "NOT_DESCENDED" if cand else ""))
            stopped, parent, parent_conf = False, None, None
            continue
        if not best or best[1] < CALL_MIN:
            out.append((None, "none", "NO_CALL" if cand else ""))
            parent, parent_conf = None, None
            continue
        gid, bs = best
        g = groups[gid]
        others = [x for x, s in at if x != gid and s >= MIXED_MIN and not related(x, gid)]
        bar = g["loo"][0] if old else g["p05"]
        conf = "high" if bs >= bar and not others else "medium" if bs >= 0.5 * g["loo"][1] else "low"
        if not old:
            if low_cov and parent_conf and RANK[conf] > RANK[parent_conf]:
                conf = parent_conf
            stopped = conf == "low" and not low_cov
        out.append((gid, conf, "MIXED" if others else ""))
        parent, parent_conf = gid, conf
    return out


def read_tsv(p):
    with open(p) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def load_groups(d):
    b = tomllib.load(open(d / "bundle" / "bundle.toml", "rb"))
    sens = defaultdict(list)
    for r in read_tsv(d / "loo.tsv"):
        sens[r["group"]].append(float(r["sensitivity"]))
    groups = {g["id"]: {"level": g["level"], "loo": g["loo"],
                        "p05": g.get("loo_p05", p05(sens[g["id"]]) if sens[g["id"]] else g["loo"][0])}
              for g in b["groups"]}
    return groups, b["levels"]


# The trained group directly above each held-out clade. Second-round members carry only a species
# label in heldout.tsv (they trained under literature labels), so the path is set per holdout.
PARENT = {"VNB": "C.neoformans", "VNII": "C.neoformans",
          "VNBI": "C.neoformans/VNB", "VNBII": "C.neoformans/VNB",
          "VNIa": "C.neoformans/VNI", "VNIb": "C.neoformans/VNI", "VNIc": "C.neoformans/VNI",
          "VNIa_C1": "C.neoformans/VNI/VNIa", "VNIa_C2": "C.neoformans/VNI/VNIa",
          "VNIa_C3": "C.neoformans/VNI/VNIa"}


def truth(tag, r, level):
    path = PARENT.get(tag, r["species"]).split("/")  # other1/other2: a deep clade of a species
    i = LEVELS.index(level)
    return "/".join(path[: i + 1]) if i < len(path) else None


def held_level(tag):
    if tag in ("VNII", "VNB"):
        return "lineage"
    if tag.startswith("VNIa_C") or tag.startswith("other"):
        return None
    return "sublineage"


def main(out, check=False):
    out = Path(out)
    train = {r["sample"] for r in read_tsv(out / "qc.tsv") if r["train"] == "true"}
    rows, mismatches, checked = [], 0, 0
    for d in sorted((out / "r2").iterdir()):
        tag = d.name
        if tag.startswith("exact_") or not (d / "heldout_calls.tsv").exists():
            continue
        groups, levels = load_groups(d)
        held = {r["sample"]: r for r in read_tsv(d / "heldout.tsv")}
        hl = held_level(tag)
        # other1/other2 hold out one clade in each of five species: report each species apart
        tally = defaultdict(lambda: defaultdict(Counter))
        n, n_lit, unexpl = Counter(), Counter(), Counter()
        score = defaultdict(lambda: defaultdict(list))
        for c in read_tsv(d / "heldout_calls.tsv"):
            s = c["sample"]
            if s not in train or s not in held:
                continue
            sc = {g: float(c["score:" + g]) for g in groups if c.get("score:" + g) not in (None, "", "NA")}
            low_cov = c["low_coverage"] == "true"
            calls = call(groups, levels, sc, low_cov, old=check)
            if check:
                for lv, (gid, conf, _) in zip(LEVELS, calls):
                    checked += 1
                    stored = (c[lv] if c[lv] != "NA" else None, c[lv + "_confidence"])
                    if stored != (gid, conf):
                        mismatches += 1
                        if mismatches <= 10:
                            print("MISMATCH", tag, s, lv, stored, (gid, conf))
                continue
            sub = held[s]["species"] if tag.startswith("other") else ""
            n[sub] += 1
            n_lit[sub] += held[s]["how"].startswith("literature") or tag.startswith("other")
            unexpl[sub] += "UNEXPLAINED" in c["flags"]
            for lv, (gid, conf, flag) in zip(LEVELS, calls):
                t = truth(tag, held[s], lv)
                if hl and LEVELS.index(lv) >= LEVELS.index(hl):
                    if LEVELS.index(lv) > LEVELS.index(hl):
                        continue  # below a held-out level there is nothing to call
                    key = "right" if gid is None else f"wrong:{gid.rsplit('/', 1)[-1]}:{conf}"
                elif t is None:
                    continue
                else:
                    if c.get("score:" + t) not in (None, "", "NA"):
                        score[sub][lv].append(float(c["score:" + t]))
                    key = f"right:{conf}" if gid == t else "no_call" if gid is None else f"wrong:{gid.rsplit('/', 1)[-1]}:{conf}"
                tally[sub][lv][key] += 1
        if check:
            continue
        for sub in sorted(tally):
            for lv in LEVELS:
                t, sc_ = tally[sub][lv], score[sub][lv]
                if t:
                    rows.append(dict(holdout=tag, species=sub, level=lv,
                                     role="held out" if hl == lv else "trained",
                                     n=n[sub], n_lit=n_lit[sub], unexplained=unexpl[sub],
                                     score_med=f"{median(sc_):.3f}" if sc_ else "",
                                     **{k: t[k] for k in sorted(t)}))
    if check:
        print(f"{checked} stored level calls, {mismatches} not reproduced")
        sys.exit(1 if mismatches else 0)
    keys = ["holdout", "species", "level", "role", "n", "n_lit", "unexplained", "score_med"]
    extra = sorted({k for r in rows for k in r} - set(keys))
    with open(out / "r2" / "table2_long.tsv", "w") as f:
        w = csv.DictWriter(f, keys + extra, delimiter="\t", restval=0)
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print("\t".join(f"{k}={r[k]}" for k in r))


if __name__ == "__main__":
    main(sys.argv[1], "--check" in sys.argv[2:])
