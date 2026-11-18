"""Coverage calibration summary: one row per run x depth x level, from census/calib/out/<run>.tsv.

Coverage is base coverage of the kept reads (spots kept x read length, both mates / the "*" solid
reference, 18.24 Mb). Confidence is recomputed with the LOW_COVERAGE parent cap (3e2770d): the runs
were classified just before it, and the cap is a function of the row alone, so this matches what
the current classify reports.

Writes census/calib/summary.tsv and prints the numbers quoted in the paper.
Usage: python census/calib/summarize.py
"""
import glob
import os

import pandas as pd

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../.."))
GENOME = 18238061.0  # census bundle solid_reference "*"
LEVELS = [("species", 5), ("lineage", 8), ("sublineage", 11)]  # call column; score, confidence follow
RANK = {"low": 0, "medium": 1, "high": 2}

runs = pd.read_csv("census/calib_runs.tsv", sep="\t", dtype=str).fillna("").set_index("sample")
rows = []
for f in sorted(glob.glob("census/calib/out/*.tsv")):
    run = os.path.basename(f)[:-4]
    t = runs.loc[run]
    want = {"species": t.species,
            "lineage": f"{t.species}/{t.lineage}" if t.lineage else "",
            "sublineage": f"{t.species}/{t.lineage}/{t.sublineage}" if t.sublineage else ""}
    kept = pd.read_csv(f"census/calib/out/{run}.spots", sep="\t", header=None, index_col=0)[1]
    tab = pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False)
    for r in tab.itertuples(index=False):
        depth = int(r[0].split("@")[1])
        low = r[2] == "true"
        parent = None
        for level, c in LEVELS:
            call, score, conf = r[c], r[c + 1], r[c + 2]
            if call != "NA" and conf != "none":
                if low and parent is not None and RANK[conf] > RANK[parent]:
                    conf = parent
                parent = conf
            else:
                parent = None
            if not want[level]:
                continue
            outcome = "no_call" if call == "NA" else ("right" if call == want[level] else "wrong")
            rows.append(dict(run=run, group=t.lineage or t.species, level=level, spots=depth,
                             kept=int(kept[depth]), read_len=int(t.rl),
                             coverage=kept[depth] * int(t.rl) / GENOME, peak=int(r[1]), low_coverage=low,
                             kmer_excess=float(r[4]), score=float(score) if score != "NA" else float("nan"),
                             outcome=outcome, confidence=conf if outcome != "no_call" else "none"))
s = pd.DataFrame(rows)
s.to_csv("census/calib/summary.tsv", sep="\t", index=False, float_format="%.4f")

print(f"{s.run.nunique()} runs, {s[s.level == 'species'].shape[0]} run x depth")
for level, g in s.groupby("level", sort=False):
    n = g.outcome.value_counts()
    hi = g[g.confidence == "high"]
    not_hi = g[g.confidence != "high"]
    print(f"{level:10s} right {n.get('right', 0):3d}  wrong {n.get('wrong', 0)}  no call {n.get('no_call', 0):3d}"
          f"  | high from {hi.coverage.min():.1f}x; non-high up to {not_hi.coverage.max():.1f}x"
          f" (peak <= {not_hi.peak.max()}); every row with peak >= 5 high: {(g[g.peak >= 5].confidence == 'high').all()}")
sp = s[s.level == "species"]
print("\nspecies, by depth: no call / medium-or-low / high, coverage range")
for d, g in sp.groupby("spots"):
    v = g.confidence.value_counts()
    print(f"  {d:>8} {v.get('none', 0):2d} {v.get('medium', 0) + v.get('low', 0):2d} {v.get('high', 0):2d}"
          f"   {g.coverage.min():5.1f}-{g.coverage.max():5.1f}x")
called = s[s.outcome == "right"]
print(f"\nscore vs kmer_excess (called rows, excess < 0.95): median |diff| "
      f"{(called[called.kmer_excess < 0.95].score - called[called.kmer_excess < 0.95].kmer_excess).abs().median():.3f}")
print("\ncoverage at the first high call per run and level (median, max):")
first = s[s.confidence == "high"].groupby(["level", "run"]).coverage.min().groupby("level")
print(first.agg(["median", "max"]).round(1).to_string())
print("\nslowest runs (first high species call):")
f = s[(s.level == "species") & (s.confidence == "high")].groupby("run").agg(cov=("coverage", "min"),
                                                                              rl=("read_len", "first"))
print(f.sort_values("cov").tail(4).round(1).to_string())
