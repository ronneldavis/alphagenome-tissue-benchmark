import sys,numpy as np,pandas as pd
SP=sys.argv[1]; G=['ISLET','PANCREAS','LIVER','ADIPOSE','MUSCLE']
T=pd.read_csv(f"{SP}/t2d_tissue_z.csv"); L=pd.read_csv(f"{SP}/ldl_tissue_z.csv")
H=pd.read_csv(f"{SP}/t2d_histone_z.csv"); HL=pd.read_csv(f"{SP}/ldl_histone_z.csv")
T=T.merge(H[['variant_id']+[g+'_hz' for g in G]],on='variant_id')
L=L.merge(HL[['variant_id']+[g+'_hz' for g in G]],on='variant_id')
for D in (T,L):
    D['other_max']=D[[g for g in G if g!='ISLET']].max(axis=1)
    D['islet_spec']=D['ISLET']-D['other_max']
    D['liver_spec']=D['LIVER']-D[[g for g in G if g!='LIVER']].max(axis=1)
    D['any_max']=D[G].max(axis=1)

print("="*94); print("T2D loci that look SPECIFICALLY islet (islet z high AND above every other tissue)"); print("="*94)
s=T[(T.ISLET>0.5)].sort_values('islet_spec',ascending=False).head(10)
print(s[['rsid','l2g_gene','pip','ISLET','ISLET_hz','LIVER','ADIPOSE','MUSCLE','islet_spec']].round(2).to_string(index=False))

print("\n"+"="*94); print("T2D loci with the biggest effect in ANY tissue (are they tissue-specific at all?)"); print("="*94)
s=T.sort_values('any_max',ascending=False).head(10)
print(s[['rsid','l2g_gene','pip','ISLET','PANCREAS','LIVER','ADIPOSE','MUSCLE','islet_spec']].round(2).to_string(index=False))

print("\n"+"="*94); print("LDL loci that look SPECIFICALLY liver (the control working)"); print("="*94)
s=L[(L.LIVER>0.5)].sort_values('liver_spec',ascending=False).head(10)
print(s[['rsid','l2g_gene','pip','LIVER','LIVER_hz','ISLET','ADIPOSE','MUSCLE','liver_spec']].round(2).to_string(index=False))

print("\n"+"="*94); print("TCF7L2 - the most replicated T2D variant on earth - what does AlphaGenome see?"); print("="*94)
print(T[T.l2g_gene=='TCF7L2'][['rsid','pip','ISLET','ISLET_hz','PANCREAS','LIVER','ADIPOSE','MUSCLE','any_max']].round(2).to_string(index=False))
print("\nfor scale, background 95th pct of any_max: n/a; T2D any_max median =",round(T.any_max.median(),2),
      "| 90th pct =",round(T.any_max.quantile(.9),2))
T.to_csv(f"{SP}/t2d_final.csv",index=False); L.to_csv(f"{SP}/ldl_final.csv",index=False)
print("\ncounts: T2D loci with any_max>1:",int((T.any_max>1).sum()),"/",len(T),
      "| LDL:",int((L.any_max>1).sum()),"/",len(L))
