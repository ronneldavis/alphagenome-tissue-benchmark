import json,sys,pandas as pd
SP=sys.argv[1]
cs=json.load(open(f"{SP}/t2d_credsets.json"))
rows=[]
for c in cs:
    mem=c["locus"]["rows"] or []
    m95=[m for m in mem if m.get("is95CredibleSet")] or mem
    if not m95: continue
    top=max(m95,key=lambda m: m.get("posteriorProbability") or 0)
    v=top["variant"]
    if not v: continue
    g=(c.get("l2GPredictions") or {}).get("rows") or []
    rows.append(dict(
        signal=c["studyLocusId"],
        variant_id=v["id"], rsid=(v.get("rsIds") or [None])[0],
        chrom=v["chromosome"], pos=v["position"],
        ref=v["referenceAllele"], alt=v["alternateAllele"],
        pip=top.get("posteriorProbability"),
        cs_size=len(m95),
        csq=(v.get("mostSevereConsequence") or {}).get("label"),
        l2g_gene=g[0]["target"]["approvedSymbol"] if g else None,
        l2g_score=round(g[0]["score"],3) if g else None,
        neglog10p=-(c["pValueExponent"]) if c.get("pValueExponent") else None,
    ))
df=pd.DataFrame(rows).drop_duplicates("variant_id")
df["coding"]=df["csq"].fillna("").str.contains("missense|stop|frameshift|synonymous|splice_acceptor|splice_donor|inframe",case=False)
print("signals:",len(df))
print("\nconsequence breakdown:"); print(df["csq"].value_counts().head(12).to_string())
print("\ncoding:",int(df.coding.sum()),"| non-coding:",int((~df.coding).sum()))
print("\nfine-mapping quality: top-PIP distribution")
print(df["pip"].describe().round(3).to_string())
print("\ncredible-set size distribution"); print(df["cs_size"].describe().round(1).to_string())
df.to_csv(f"{SP}/t2d_variants.csv",index=False)
print("\n--- top 15 best-resolved non-coding signals ---")
s=df[~df.coding].sort_values("pip",ascending=False).head(15)
print(s[["rsid","chrom","pos","ref","alt","pip","cs_size","l2g_gene","csq"]].to_string(index=False))
