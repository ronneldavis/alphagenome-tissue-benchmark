import sys,json,numpy as np,pandas as pd
SP=sys.argv[1]; G=['ISLET','PANCREAS','LIVER','ADIPOSE','MUSCLE']
T=pd.read_csv(f"{SP}/t2d_final.csv"); L=pd.read_csv(f"{SP}/ldl_final.csv")
t=T.sort_values('any_max',ascending=False)
out=dict(
 t2d_loci=[dict(rsid=r.rsid,gene=r.l2g_gene,v=round(float(r.any_max),2),
                islet=round(float(r.ISLET),2)) for r in t.itertuples()],
 n_t2d=len(T), n_ldl=len(L),
 t2d_above1=int((T.any_max>1).sum()), ldl_above1=int((L.any_max>1).sum()),
 tcf7l2=[dict(rsid=r.rsid,pip=float(r.pip),islet=round(float(r.ISLET),2),
              anymax=round(float(r.any_max),2)) for r in T[T.l2g_gene=='TCF7L2'].itertuples()],
)
json.dump(out,open(f"{SP}/chart_data.json","w"))
print("t2d loci exported:",len(out['t2d_loci']),"| above1:",out['t2d_above1'],"/",out['n_t2d'])
print("top6:",[(d['gene'],d['v']) for d in out['t2d_loci'][:6]])
print("tcf7l2:",out['tcf7l2'])
print("rank of TCF7L2 rs7903146 by any_max:", int(t.reset_index().index[t.reset_index().rsid=='rs7903146'][0])+1, "of", len(t))
print("median any_max:",round(float(T.any_max.median()),2))
