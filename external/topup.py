"""Pick reserve isolates for pure lineages that have fewer than 12 good databases.

Good = a presence database exists and is not flagged low_cov in stats.tsv (no genome peak, or a peak
below 5x). Old runs can be so error-rich that 1M pairs gives no usable genome peak; this can't be
seen in the metadata, so reserves replace them. Writes topup.tsv (next reserves, in order, that have
no database yet) and prints the shortfall per lineage; exits 1 when nothing is needed.

Usage: python topup.py
"""
import os
import sys
import pandas as pd

N_PURE = 12
os.chdir(os.path.dirname(os.path.abspath(__file__)))
s = pd.read_csv("samples.tsv", sep="\t")
st = pd.read_csv("stats.tsv", sep="\t").set_index("run") if os.path.exists("stats.tsv") else pd.DataFrame()
has_db = s.run.map(lambda r: os.path.exists(f"db/{r}.kmc_suf"))
good = has_db & s.run.map(lambda r: r in st.index and st.low_cov.get(r) == 0)
tried = has_db | s.run.isin(open("logs/fetch.log").read().split() if os.path.exists("logs/fetch.log") else [])
pick = []
for lin, g in s[s.set != "hybrid"].groupby("lineage"):
    need = N_PURE - int(good[g.index].sum())
    spare = g[(g.set == "reserve") & ~tried[g.index]]
    print(f"{lin}: {int(good[g.index].sum())} good, need {max(0, need)}, {len(spare)} reserves left")
    if need > 0:
        pick.append(spare.head(need))
pick = pd.concat(pick) if pick else s.head(0)
pick.to_csv("topup.tsv", sep="\t", index=False)
sys.exit(0 if len(pick) else 1)
