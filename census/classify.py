import csv,glob,collections,json,re
rows=[r for f in glob.glob("runinfo_*.csv") for r in csv.DictReader(open(f)) if (r.get("Run") or "")[:3] in("SRR","ERR","DRR")]
bs=json.load(open("biosamples.json")); P=json.load(open("projects.json"))["projects"]
EXP=r'adaptor|mutant|evolution|evolv|H99|KN99|JEC2|XL280|deletion|treated|treatment|hypermut|RNAi|centromere|profiling|library prep|reference|calcineurin|ChIP|uniparental|transposon|transposition|thermotoleran|CO2|tolerance|revertant|resistan|selective|amoeba|macrophage|MHC|statin|mechanism|multiomics|knockout|aneuploid|transcription factor|drug'
SERIAL=r'serial|reoccurrence|recurren|relapse|over the course'
LAB=r'^(H99|KN99|JEC2[01]|XL280|B-?3501|CM\d|R265|WM276)|H99|KN99|_mut|Δ|delta|::|adaptor|evolved|passage|clone ?\d|ancestor'
def species(org):
    o=org.lower()
    for s in ["neoformans","deneoformans","deuterogattii","bacillisporus","tetragattii","decagattii","hyracis","gattii","brasiliensis","schenckii","globosa","luriei"]:
        if s in o: return s
    return "other_"+org.split()[0]
recs={}
for r in rows:
    b=r["BioSample"]
    if b in recs: continue
    a=bs.get(b,{}); p=r["BioProject"]; t=P.get(p,{}).get("title","")
    blob=" ".join([a.get("strain",""),a.get("isolate",""),a.get("sample_name",""),a.get("_title",""),a.get("_org",""),a.get("genotype",""),a.get("karyotype",""),a.get("sample_type","")])
    org=a.get("_org") or r["ScientificName"]
    sp=species(org)
    # deneoformans often filed under C. neoformans var. neoformans
    if re.search(r'var\.? neoformans|deneoformans|serotype D',org+blob,re.I) and sp=="neoformans": sp="deneoformans"
    if a.get("serotype")=="D": sp="deneoformans"
    if re.search(r'AD hybrid|serotype AD|\bAD\b',blob) or a.get("serotype")=="AD" or "hybrid" in t.lower(): cls="hybrid/deneo-study"
    elif re.search(SERIAL,t,re.I): cls="serial"
    elif re.search(EXP,t,re.I) or re.search(LAB,blob,re.I) or re.search(r'aneuploid|mutant|null',blob,re.I): cls="lab/experimental"
    else: cls="natural"
    strain=(a.get("strain") or a.get("isolate") or a.get("sample_name") or r["SampleName"]).strip()
    recs[b]=dict(bs=b,proj=p,sp=sp,cls=cls,strain=strain,country=a.get("geo_loc_name","").split(":")[0],
                 source=a.get("isolation_source",""),host=a.get("host",""),center=r["CenterName"],title=t)
json.dump(recs,open("classified.json","w"))
# natural: dedupe by normalized strain name
def norm(s): return re.sub(r'[^a-z0-9]','',s.lower())
nat=[x for x in recs.values() if x["cls"]=="natural"]
seen={}; dup=0
for x in nat:
    k=(x["sp"],norm(x["strain"]))
    if k[1] in("","missing","notapplicable") : k=(x["sp"],x["bs"])
    if k in seen: dup+=1; x["dup_of"]=seen[k]
    else: seen[k]=x["bs"]
print("class x species (all biosamples):")
tab=collections.Counter((x["sp"],x["cls"]) for x in recs.values())
sps=sorted({s for s,_ in tab}); cl=["natural","serial","lab/experimental","hybrid/deneo-study"]
print(f'{"species":18}'+"".join(f'{c[:12]:>14}' for c in cl)+"   natural_dedup")
for s in sps:
    nd=sum(1 for x in nat if x["sp"]==s and "dup_of" not in x)
    print(f'{s:18}'+"".join(f'{tab[(s,c)]:14}' for c in cl)+f'{nd:14}')
print("dup strain names collapsed:",dup)
json.dump(recs,open("classified.json","w"))
