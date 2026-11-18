#!/usr/bin/env python3
"""Add loo_p05 (5th-percentile leave-one-out sensitivity) to a bundle built before cladeprint
wrote it, from the run's loo.tsv. classify uses it as the bar for a high-confidence call.

Same nearest-rank quantile as pipeline.rs (index round((n - 1) * 0.05) of the sorted values).
loo.tsv rounds to 4 decimals, so each group's min and median are checked against bundle.toml
(to 1e-4) before anything is written: a loo.tsv from another run fails the check.

usage: add_loo_p05.py BUNDLE_DIR LOO_TSV
"""
import csv
import re
import sys
import tomllib
from collections import defaultdict
from pathlib import Path


def main(bdir, loo_tsv):
    path = Path(bdir) / "bundle.toml"
    text = path.read_text()
    groups = tomllib.loads(text)["groups"]
    sens = defaultdict(list)
    with open(loo_tsv) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            sens[r["group"]].append(float(r["sensitivity"]))
    p05 = {}
    for g in groups:
        v = sorted(sens.get(g["id"], []))
        if not v:
            sys.exit(f"{g['id']}: no rows in {loo_tsv}")
        n = len(v)
        med = v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2
        if abs(v[0] - g["loo"][0]) > 1e-4 or abs(med - g["loo"][1]) > 1e-4:
            sys.exit(f"{g['id']}: loo.tsv min/median {v[0]}/{med} != bundle {g['loo'][:2]}")
        p05[g["id"]] = v[round((n - 1) * 0.05)]
    # insert after each group's loo line; groups are matched by their id line
    out, gid = [], None
    for line in text.splitlines(keepends=True):
        if line.startswith("loo_p05"):
            continue  # rerun: replace
        out.append(line)
        m = re.match(r'id = "(.*)"', line)
        if m:
            gid = m.group(1)
        if line.startswith("loo = ") and gid is not None:
            out.append(f"loo_p05 = {p05[gid]}\n")
            gid = None
    path.write_text("".join(out))
    for g in groups:
        print(f"{g['id']}\tn={len(sens[g['id']])}\tmin={g['loo'][0]:.4f}\tp05={p05[g['id']]:.4f}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
