import csv,glob,time,urllib.request,urllib.parse,xml.etree.ElementTree as ET,json
E="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
acc=set()
for f in glob.glob("runinfo_*.csv"):
    for r in csv.DictReader(open(f)):
        if (r.get("BioSample") or "").startswith("SAM"): acc.add(r["BioSample"])
acc=sorted(acc); print(len(acc),"biosamples")
out={}
def post(url,data):
    for i in range(5):
        try: return urllib.request.urlopen(url,urllib.parse.urlencode(data).encode(),timeout=120).read()
        except Exception as e: print("retry",e); time.sleep(3)
for i in range(0,len(acc),200):
    chunk=acc[i:i+200]
    x=post(f"{E}/esearch.fcgi",{"db":"biosample","term":" OR ".join(c+"[accn]" for c in chunk),"retmax":500,"usehistory":"y"})
    r=ET.fromstring(x); we=r.findtext("WebEnv"); qk=r.findtext("QueryKey")
    x=post(f"{E}/efetch.fcgi",{"db":"biosample","query_key":qk,"WebEnv":we,"retmax":500})
    for bs in ET.fromstring(x).findall("BioSample"):
        a={at.get("harmonized_name") or at.get("attribute_name"):(at.text or "") for at in bs.iter("Attribute")}
        a["_title"]=bs.findtext("Description/Title") or ""
        a["_org"]=bs.findtext("Description/Organism/OrganismName") or ""
        out[bs.get("accession")]=a
    print(i,len(out),flush=True); time.sleep(0.5)
json.dump(out,open("biosamples.json","w"))
