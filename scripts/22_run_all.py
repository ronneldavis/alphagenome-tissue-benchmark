import os,sys,json,time,urllib.request,numpy as np,pandas as pd
from alphagenome.models import dna_client, variant_scorers
from alphagenome.data import genome
SP=sys.argv[1]; CAP=150
URL="https://api.platform.opentargets.org/api/v4/graphql"
def gq(s,tries=4):
    for t in range(tries):
        try:
            r=urllib.request.Request(URL,data=json.dumps({"query":s}).encode(),headers={"Content-Type":"application/json"})
            return json.load(urllib.request.urlopen(r,timeout=180))
        except Exception:
            if t==tries-1: raise
            time.sleep(3)
Q="""{ credibleSets(studyIds:["%s"], page:{index:%d,size:25}){ count rows{
 studyLocusId pValueExponent
 locus(page:{index:0,size:200}){ rows{ posteriorProbability is95CredibleSet
   variant{ id chromosome position referenceAllele alternateAllele rsIds mostSevereConsequence{ label } } } }
 l2GPredictions(page:{index:0,size:1}){ rows{ score target{ approvedSymbol } } } } } }"""

panel=json.load(open(f"{SP}/trait_panel.json"))
panel.append(dict(trait="type 2 diabetes (replication)",expected="ISLET",study="GCST006867",dname="T2D Xue 2018"))
m=dna_client.create(os.environ['ALPHAGENOME_KEY'])
W=dna_client.SUPPORTED_SEQUENCE_LENGTHS['SEQUENCE_LENGTH_1MB']
SC=[variant_scorers.RECOMMENDED_VARIANT_SCORERS[k] for k in ('RNA_SEQ','CHIP_HISTONE')]
CODING='missense|stop|frameshift|synonymous|splice_acceptor|splice_donor|inframe'

for tr in panel:
    tag=tr['study']
    if os.path.exists(f"{SP}/T_{tag}_scores.npz"): print("skip",tag,flush=True); continue
    rows=[];i=0
    while True:
        d=gq(Q%(tag,i))["data"]["credibleSets"]; rows+=d["rows"]; i+=1
        if len(rows)>=d["count"] or not d["rows"] or len(rows)>=400: break
    recs=[]
    for c in rows:
        mem=[x for x in (c["locus"]["rows"] or []) if x.get("is95CredibleSet")] or (c["locus"]["rows"] or [])
        if not mem: continue
        top=max(mem,key=lambda x:x.get("posteriorProbability") or 0); v=top["variant"]
        if not v: continue
        g=(c.get("l2GPredictions") or {}).get("rows") or []
        recs.append(dict(variant_id=v["id"],rsid=(v.get("rsIds") or [None])[0],chrom=v["chromosome"],
            pos=v["position"],ref=v["referenceAllele"],alt=v["alternateAllele"],
            pip=top.get("posteriorProbability"),cs_size=len(mem),
            csq=(v.get("mostSevereConsequence") or {}).get("label"),
            l2g_gene=g[0]["target"]["approvedSymbol"] if g else None,
            negp=-(c.get("pValueExponent") or 0)))
    df=pd.DataFrame(recs).drop_duplicates("variant_id")
    df=df[~df["csq"].fillna("").str.contains(CODING,case=False)]
    df=df.sort_values("negp",ascending=False).head(CAP).reset_index(drop=True)
    df.to_csv(f"{SP}/T_{tag}_variants.csv",index=False)
    rr,hh,meta=[],[],[];t0=time.time()
    for _,r in df.iterrows():
        v=genome.Variant(chromosome=f"chr{r.chrom}",position=int(r.pos),
                         reference_bases=str(r.ref),alternate_bases=str(r.alt))
        out=None
        for a in range(3):
            try: out=m.score_variant(interval=v.reference_interval.resize(W),variant=v,variant_scorers=SC); break
            except Exception: time.sleep(2*(a+1))
        if out is None: continue
        X=np.asarray(out[0].X); genes=out[0].obs['gene_name'].astype(str).tolist()
        gi=genes.index(r.l2g_gene) if r.l2g_gene in genes else (int(np.nanargmax(np.nanmean(np.abs(X),axis=1))) if X.size else 0)
        rr.append(X[gi]); hh.append(np.asarray(out[1].X)[0])
        meta.append(dict(variant_id=r.variant_id,rsid=r.rsid,pip=r.pip,l2g_gene=r.l2g_gene,cs_size=r.cs_size))
    np.savez_compressed(f"{SP}/T_{tag}_scores.npz",rna=np.array(rr),hist=np.array(hh))
    pd.DataFrame(meta).to_csv(f"{SP}/T_{tag}_meta.csv",index=False)
    print(f"{tr['trait'][:32]:34s} {tag:<16s} n={len(rr):<4d} {time.time()-t0:.0f}s",flush=True)
json.dump(panel,open(f"{SP}/trait_panel_full.json","w"))
print("ALL DONE",flush=True)
