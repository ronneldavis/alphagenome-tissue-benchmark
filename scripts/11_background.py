import sys,json,random,urllib.request,urllib.error,time
import pandas as pd
SP=sys.argv[1]; random.seed(7)
tv=pd.concat([pd.read_csv(f"{SP}/t2d_variants.csv"),pd.read_csv(f"{SP}/ldl_variants.csv")])
tv=tv[~tv.coding]
# place 2 background SNVs per real signal, 50-400kb away, same chromosome/neighbourhood
cand=[]
for _,r in tv.iterrows():
    for _ in range(2):
        off=random.randint(50_000,400_000)*random.choice([-1,1])
        p=int(r.pos)+off
        if p>1_000_000: cand.append((str(r.chrom),p))
random.shuffle(cand); cand=cand[:320]
print("candidate background positions:",len(cand))

# batch-fetch reference base from Ensembl GRCh38
URL="https://rest.ensembl.org/sequence/region/human"
out=[]
for i in range(0,len(cand),50):
    chunk=cand[i:i+50]
    body={"regions":[f"{c}:{p}..{p}:1" for c,p in chunk]}
    req=urllib.request.Request(URL,data=json.dumps(body).encode(),
        headers={"Content-Type":"application/json","Accept":"application/json"})
    for a in range(4):
        try:
            res=json.load(urllib.request.urlopen(req,timeout=120)); break
        except Exception as e:
            if a==3: print("ensembl fail",str(e)[:80]); res=[]; break
            time.sleep(3)
    for (c,p),s in zip(chunk,res):
        b=(s.get("seq") or "").upper()
        if b in "ACGT": out.append(dict(chrom=c,pos=p,ref=b))
    print(f"  fetched {len(out)}",flush=True)
    time.sleep(0.3)
TS={'A':'G','G':'A','C':'T','T':'C'}
bg=pd.DataFrame(out); bg['alt']=bg['ref'].map(TS)
bg['variant_id']=bg.chrom.astype(str)+"_"+bg.pos.astype(str)+"_"+bg.ref+"_"+bg.alt
bg['rsid']=None; bg['pip']=None; bg['cs_size']=None; bg['l2g_gene']=None; bg['coding']=False
bg.to_csv(f"{SP}/bg_variants.csv",index=False)
print("background variants:",len(bg)); print(bg.head(3).to_string())
