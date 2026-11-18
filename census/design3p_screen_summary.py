import pandas as pd
o="/home/ubuntu/cladeprint/census/out/"
s=pd.read_csv(o+"design3p/sites_top.tsv",sep="\t",dtype={"mismatches":str})
s["mm"]=pd.to_numeric(s.mismatches,errors="coerce"); s["t"]=pd.to_numeric(s.mm_3prime,errors="coerce")
s["binds"]=s.mm.le(3)&s.t.eq(0)
ol=pd.read_csv(o+"design3p/oligos_top.tsv",sep="\t"); rep={}
for r in ol.itertuples(): rep.setdefault(r.sequence,r.oligo)
p=pd.read_csv(o+"design3p/pairs_top.tsv",sep="\t")
lab=pd.read_csv(o+"subclade/samples.tsv",sep="\t",dtype=str).set_index("sample")
qc=pd.read_csv(o+"qc.tsv",sep="\t",dtype=str).set_index("sample")
lab=lab[qc.reindex(lab.index).train.values=="true"]
L=["species","lineage","sublineage","subclade"]
def gp(lv):
    c=L[:L.index(lv)+1]
    return lab[c].apply(lambda r: None if any(pd.isna(x) or x=="" for x in r) else "/".join(r),axis=1)
mem={}
for lv in L:
    pa=gp(lv)
    for g in pa.dropna().unique(): mem[g]=set(pa.index[pa==g])
by=s.set_index(["oligo","sample"]).binds
seqof=dict(zip(ol.oligo,ol.sequence))
rows=[]
for r in p.itertuples():
    g=r.group.replace("_","/")
    g=[k for k in mem if k.replace("/","_")==r.group]
    if not g: continue
    g=g[0]; m=mem[g]; non=[x for x in lab.index if x not in m]
    L2=rep.get(r.left); R2=rep.get(r.right)
    if L2 not in by.index.get_level_values(0) or R2 not in by.index.get_level_values(0): continue
    ok=by.loc[L2].reindex(lab.index).fillna(False)&by.loc[R2].reindex(lab.index).fillna(False)
    rows.append(dict(group=g,pair=r.pair,product_bp=r.product_bp,
        member=round(ok.reindex(list(m)).mean(),4),nonmember=round(ok.reindex(non).mean(),4),
        left=r.left,right=r.right))
d=pd.DataFrame(rows); d["margin"]=d.member-d.nonmember
d.to_csv(o+"design3p/pairs_screened.tsv",sep="\t",index=False)
best=d.sort_values("margin",ascending=False).groupby("group").head(1)
print(best[["group","product_bp","member","nonmember","margin"]].sort_values("margin",ascending=False).to_string(index=False))
print("\npairs with member>=0.95 and nonmember<=0.01, per group:")
print(d[(d.member>=0.95)&(d.nonmember<=0.01)].groupby("group").size().to_string())
