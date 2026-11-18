import csv,glob,json,time,collections,urllib.request,urllib.parse,xml.etree.ElementTree as ET
E="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
def get(p,**q):
    for i in range(5):
        try: return urllib.request.urlopen(f"{E}/{p}?"+urllib.parse.urlencode(q),timeout=60).read()
        except Exception as e: time.sleep(3)
bs=json.load(open("biosamples.json"))
proj=collections.defaultdict(lambda:{"bs":set(),"runs":0,"org":collections.Counter(),"pmid":set(),"bases":0})
for f in glob.glob("runinfo_*.csv"):
    for r in csv.DictReader(open(f)):
        if not (r.get("Run") or "")[:3] in("SRR","ERR","DRR"): continue
        p=proj[r["BioProject"]]; p["bs"].add(r["BioSample"]); p["runs"]+=1; p["org"][r["ScientificName"]]+=1
        if r.get("Study_Pubmed_id") not in (None,"","0"): p["pmid"].add(r["Study_Pubmed_id"])
        try: p["bases"]+=int(r["bases"] or 0)
        except: pass
big=[k for k,v in proj.items() if len(v["bs"])>=5]
print(len(proj),"projects;",len(big),"with >=5 biosamples")
res={}
for k in big:
    x=get("esearch.fcgi",db="bioproject",term=k); uid=ET.fromstring(x).findtext("IdList/Id")
    title=""; pm=set(proj[k]["pmid"])
    if uid:
        j=json.loads(get("esummary.fcgi",db="bioproject",id=uid,retmode="json"))["result"][uid]
        title=j.get("project_title","")
        l=ET.fromstring(get("elink.fcgi",dbfrom="bioproject",db="pubmed",id=uid))
        pm|={i.text for i in l.iter("Id")}-{uid}
    v=proj[k]
    res[k]=dict(n_bs=len(v["bs"]),runs=v["runs"],org=v["org"].most_common(3),title=title,pmids=sorted(pm),
               gb_per_sample=round(v["bases"]/1e9/max(1,len(v["bs"])),2))
    time.sleep(0.35)
allpm=sorted({p for r in res.values() for p in r["pmids"]})
cites={}
for i in range(0,len(allpm),150):
    j=json.loads(get("esummary.fcgi",db="pubmed",id=",".join(allpm[i:i+150]),retmode="json"))["result"]
    for u in j.get("uids",[]): cites[u]=f'{(j[u].get("sortfirstauthor") or "")} {j[u].get("pubdate","")[:4]} {j[u].get("source","")}: {j[u].get("title","")[:90]}'
json.dump({"projects":res,"pubmed":cites},open("projects.json","w"),indent=1)
