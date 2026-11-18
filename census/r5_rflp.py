"""In-silico URA5-RFLP and PLB1-RFLP typing of the census.

Neither assay publishes fragment sizes (types are read off gels against standard strains), so the
amplicon is reconstructed from each genome's unitigs where both primers lie on one unitig, and
digested computationally. Patterns are then compared with the k-mer group assignment.
URA5-RFLP: HhaI (GCG^C) + Sau96I (G^GNCC) double digest (Meyer et al. 2003).
PLB1-RFLP: AvaI (C^YCGRG) and HindIII (A^AGCTT) separately (Latouche et al. 2003).
"""
import re, sys
import pandas as pd
O="/home/ubuntu/cladeprint/census/out/"
def read_fasta(p):
    n,s=None,[]
    for l in open(p):
        if l[0]==">":
            if n: yield n,"".join(s)
            n,s=l[1:].split()[0],[]
        else: s.append(l.strip())
    if n: yield n,"".join(s)
ENZ={"HhaI":("GCGC",3),"Sau96I":("GG[ACGT]CC",1),"AvaI":("C[CT]CG[AG]G",1),"HindIII":("AAGCTT",1)}
def digest(seq,enzymes):
    cuts=set()
    for e in enzymes:
        pat,off=ENZ[e]
        for mt in re.finditer(f"(?={pat})",seq): cuts.add(mt.start()+off)
    pos=[0]+sorted(cuts)+[len(seq)]
    return sorted([pos[i+1]-pos[i] for i in range(len(pos)-1)],reverse=True)
sites=pd.read_csv(O+"r5_sites.tsv",sep="\t",dtype={"mismatches":str})
lab=pd.read_csv(O+"subclade/samples.tsv",sep="\t",dtype=str).set_index("sample")
out=[]
for assay,(fo,ro,enz) in {"URA5":("A098|URA5","A099|SJ101",["HhaI","Sau96I"]),
                          "PLB1":("A114|IDPLB1","A115|IDPLB1R",["AvaI"])}.items():
    A=sites[sites.oligo==fo].set_index("sample"); B=sites[sites.oligo==ro].set_index("sample")
    if A.empty or B.empty: print(f"{assay}: oligos missing",file=sys.stderr); continue
    j=A.join(B,rsuffix="_r",how="inner")
    j=j[(j.mismatches!="none")&(j.mismatches_r!="none")&(j.contig==j.contig_r)]
    n=0
    for smp,r in j.iterrows():
        try: g=dict(read_fasta(f"{O}bench/genomes/{smp}.fa"))
        except FileNotFoundError: continue
        c=g.get(r.contig)
        if c is None: continue
        a,b=sorted([int(r.pos),int(r.pos_r)])
        amp=c[a:b+len(r.oligo)]
        if not (400<len(amp)<2500): continue
        out.append(dict(assay=assay,sample=smp,amp_len=len(amp),
                        pattern="-".join(str(x) for x in digest(amp,enz))))
        n+=1
    print(f"{assay}: {n} amplicons digested",file=sys.stderr)
d=pd.DataFrame(out)
d["species"]=d["sample"].map(lab.species); d["lineage"]=d["sample"].map(lab.lineage)
d["type"]=d.species.where(d.species!="C.neoformans",d.species+"/"+d.lineage.fillna("?"))
d.to_csv(O+"r5_rflp.tsv",sep="\t",index=False)
for assay,g in d.groupby("assay"):
    print(f"\n== {assay}: {g.pattern.nunique()} distinct patterns in {len(g)} genomes")
    ct=pd.crosstab(g.type,g.pattern)
    keep=ct.sum(0).sort_values(ascending=False).head(8).index
    print(ct[keep].to_string())
