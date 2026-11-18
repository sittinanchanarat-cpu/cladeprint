"""Pick pilot isolates: lineage-labelled anchors + a cohort-stratified random sample."""
import csv,json,random,collections
random.seed(20261011)
R=json.load(open("../census/classified.json")); bs=json.load(open("../census/biosamples.json"))
runs=collections.defaultdict(list)
for f in ["../census/runinfo_5207.csv"]:
    for r in csv.DictReader(open(f)):
        if (r.get("Run") or "")[:3] in("SRR","ERR","DRR") and r["LibraryLayout"]=="PAIRED" and r["Platform"] in("ILLUMINA","BGISEQ","DNBSEQ"):
            runs[r["BioSample"]].append(r)
ANCH={"SAMN18057101":"VNI","SAMN15273642":"VNI","SAMN18057102":"VNBI","SAMN18057103":"VNBI","SAMN39931093":"VNBI",
      "SAMN23424732":"VNII","SAMN08284936":"VNII","SAMN49712038":"VNII","SAMN49712034":"VNII","SAMN49712036":"VNII",
      "SAMN49712039":"VNIa","SAMN49712011":"VNIa"}
vnib=[b for b,a in bs.items() if a.get("subtype")=="VNIb" and b in runs][:2]
ANCH.update({b:"VNIb" for b in vnib})
STRAT={"PRJEB5282":4,"PRJEB27222":4,"PRJNA269039":4,"PRJEB14735":4,"PRJNA764746":3,"PRJEB76062":3,"PRJNA1182044":3,
       "PRJNA1006382":2,"PRJNA549026":2,"PRJEB13814":2,"PRJNA702892":2,"PRJEB19514":3,"_BROAD":8}
P=json.load(open("../census/projects.json"))["projects"]
pool=collections.defaultdict(list)
for b,x in R.items():
    if x["sp"]=="neoformans" and x["cls"]=="natural" and x.get("indep") and b in runs and b not in ANCH:
        pool[x["proj"] if x["proj"] in P else ("_BROAD" if x["center"]=="BI" else "_other")].append(b)
sel=[(b,"anchor",l) for b,l in ANCH.items() if b in runs]
for p,n in STRAT.items(): sel+= [(b,"pilot","") for b in random.sample(sorted(pool[p]),min(n,len(pool[p])))]
with open("pilot_samples.tsv","w") as o:
    o.write("biosample\trun\trole\tlabel\tbioproject\tstrain\tcountry\tspots\tavgLength\n")
    for b,role,l in sel:
        r=max(runs[b],key=lambda r:int(r["spots"] or 0))
        o.write(f'{b}\t{r["Run"]}\t{role}\t{l}\t{r["BioProject"]}\t{R[b]["strain"]}\t{R[b]["country"]}\t{r["spots"]}\t{r["avgLength"]}\n')
print(len(sel),"samples;",sum(1 for s in sel if s[1]=="anchor"),"anchors")
