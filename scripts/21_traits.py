import json,urllib.request,urllib.error,sys,time
URL="https://api.platform.opentargets.org/api/v4/graphql"
def q(s,tries=3):
    for t in range(tries):
        try:
            r=urllib.request.Request(URL,data=json.dumps({"query":s}).encode(),headers={"Content-Type":"application/json"})
            return json.load(urllib.request.urlopen(r,timeout=120))
        except Exception as e:
            if t==tries-1: raise
            time.sleep(2)
WANT=[("type 2 diabetes","ISLET"),("LDL cholesterol","LIVER"),("HDL cholesterol","LIVER"),
 ("triglyceride measurement","LIVER"),("Crohn's disease","IMMUNE"),("rheumatoid arthritis","IMMUNE"),
 ("type 1 diabetes","IMMUNE"),("atrial fibrillation","HEART"),("coronary artery disease","ARTERY"),
 ("asthma","IMMUNE"),("body mass index","BRAIN"),("ulcerative colitis","IMMUNE"),
 ("systolic blood pressure","ARTERY"),("alkaline phosphatase measurement","LIVER"),
 ("atopic eczema","SKIN"),("schizophrenia","BRAIN"),("FEV/FEC ratio","LUNG"),
 ("heart rate","HEART"),("alanine aminotransferase measurement","LIVER"),("platelet count","IMMUNE")]
res=[]
for term,exp in WANT:
    try:
        h=q('{ search(queryString:"%s", entityNames:["disease"], page:{index:0,size:3}){ hits{ id name } } }'%term)["data"]["search"]["hits"]
    except Exception as e:
        print("search fail",term); continue
    if not h: print("no hit:",term); continue
    did,dname=h[0]["id"],h[0]["name"]
    try:
        st=q('{ studies(diseaseIds:["%s"], page:{index:0,size:120}){ rows{ id studyType nSamples } } }'%did)["data"]["studies"]["rows"]
    except Exception: continue
    st=[x for x in st if x["studyType"]=="gwas" and x.get("nSamples")]
    st.sort(key=lambda x:x["nSamples"],reverse=True)
    best=None
    for x in st[:8]:
        try: c=q('{ credibleSets(studyIds:["%s"], page:{index:0,size:1}){ count } }'%x["id"])["data"]["credibleSets"]["count"]
        except Exception: c=0
        if best is None or c>best[1]: best=(x["id"],c,x["nSamples"])
        if c>=150: break
    if best and best[1]>=40:
        res.append(dict(trait=term,expected=exp,disease=did,dname=dname,study=best[0],n_cs=best[1],n=best[2]))
        print(f"{term:38s} {exp:9s} {best[0]:<20s} cs={best[1]:<4d} n={best[2]}",flush=True)
    else:
        print(f"{term:38s} SKIP (best cs={best[1] if best else 0})",flush=True)
json.dump(res,open(f"{sys.argv[1]}/trait_panel.json","w"))
print("\nusable traits:",len(res))
