"""Why the census URA5 primer audit (0.948 of C. neoformans) is lower than the assemblies (0.986).

Classifies each census C. neoformans genome's URA5/SJ101 binding state from r5_sites.tsv (best site
in the genome's unitigs, strict 3' rule as in r5_summary.py): both bind / a primer has no site within
3 mismatches ("absent") / a site with a 3'-tail mismatch ("tail"), then crosses it with k-mer
coverage (qc.tsv peak) and lineage. Writes tmp/ura5_census_state.tsv.
Usage: python census/r5_ura5_check.py
"""
import pandas as pd
sites = pd.read_csv("census/out/r5_sites.tsv", sep="\t", dtype=str)
s = pd.read_csv("census/out/prop/samples.tsv", sep="\t", dtype=str).set_index("sample")
qc = pd.read_csv("census/out/qc.tsv", sep="\t", dtype=str).set_index("sample")
neo = [x for x in s.index[s.species=="C.neoformans"] if qc.train.get(x)=="true"]
u = sites[sites.oligo.isin(["A098|URA5","A099|SJ101"])].copy()
u["mm"]=pd.to_numeric(u.mismatches,errors="coerce"); u["t"]=pd.to_numeric(u.mm_3prime,errors="coerce")
def state(r):
    if pd.isna(r.mm): return "absent"
    return "bind" if r.t==0 else f"tail{int(r.t)}"
u["st"]=u.apply(state,axis=1)
w = u.pivot(index="sample",columns="oligo",values="st").reindex(neo).fillna("absent")
w.columns=["URA5","SJ101"]
w["ok"]=(w.URA5=="bind")&(w.SJ101=="bind")
print("n",len(w),"ok",w.ok.mean().round(4))
f=w[~w.ok].copy()
print(f.groupby(["URA5","SJ101"]).size().sort_values(ascending=False))
q=qc.reindex(w.index); w["peak"]=pd.to_numeric(q.peak); w["lowcov"]=q.low_coverage
w["err_ratio"]=pd.to_numeric(q.error_kmers)/pd.to_numeric(q.solid_kmers)
print(w.groupby("ok")[["peak","err_ratio"]].median())
print(pd.crosstab(pd.cut(w.peak,[0,4,6,8,12,1000]),w.ok,normalize="index").round(3), pd.crosstab(pd.cut(w.peak,[0,4,6,8,12,1000]),w.ok))
print(pd.crosstab(s.reindex(w.index).lineage.fillna("?"), w.ok))
w.to_csv("tmp/ura5_census_state.tsv",sep="\t")
