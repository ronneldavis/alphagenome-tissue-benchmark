import sys,numpy as np,pandas as pd
from scipy.stats import mannwhitneyu, false_discovery_control
SP=sys.argv[1]
MARKS=['H3K27ac','H3K36me3','H3K4me1','H3K4me3','H3K9ac']
h=pd.read_csv(f"{SP}/tracks_chip_histone.csv")
GR=['ISLET','LIVER','ADIPOSE','MUSCLE','BRAIN','IMMUNE','HEART','ARTERY','LUNG','INTESTINE','SKIN']
ok=h['histone_mark'].isin(MARKS).values
IDX={g:((h['grp']==g).values&ok) for g in GR}
BG=np.abs(np.load(f"{SP}/bg_scores.npz")['hist'])
NR={g:(BG[:,IDX[g]].mean(1).mean(),BG[:,IDX[g]].mean(1).std()) for g in GR}

R=pd.read_csv(f"{SP}/bench_summary.csv")
q=false_discovery_control(R['p_exp'].values,method='bh')
R['q_exp']=q
print("=== expected tissue enriched vs background ===")
print(f"p<0.05: {(R.p_exp<0.05).sum()}/20 | BH q<0.05: {(R.q_exp<0.05).sum()}/20 | Bonferroni: {(R.p_exp<0.05/20).sum()}/20")
print(R.sort_values('p_exp')[['trait','expected','rank','exp_z','p_exp','q_exp']].to_string(index=False,float_format=lambda v:f"{v:.3g}"))
R.to_csv(f"{SP}/bench_summary.csv",index=False)

print("\n=== T2D: two datasets compared ===")
for sid,lab in (("GCST009379","Mahajan 2018 (n=898k)"),("GCST006867","Xue 2018 (n=659k)")):
    H=np.abs(np.load(f"{SP}/T_{sid}_scores.npz")['hist']) if sid=="GCST009379" else np.abs(np.load(f"{SP}/t2d_scores.npz")['hist'])
    z={g:(H[:,IDX[g]].mean(1)-NR[g][0])/NR[g][1] for g in GR}
    o=sorted(GR,key=lambda g:-z[g].mean())
    print(f"\n{lab}  n={H.shape[0]}")
    print("  " + " | ".join(f"{g}:{z[g].mean():+.2f}" for g in o[:5]))
    print(f"  ISLET mean {z['ISLET'].mean():+.3f} median {np.median(z['ISLET']):+.3f} "
          f"frac>1 {(z['ISLET']>1).mean()*100:.0f}%  p={mannwhitneyu(z['ISLET'],(BG[:,IDX['ISLET']].mean(1)-NR['ISLET'][0])/NR['ISLET'][1],alternative='greater').pvalue:.3f}")

vm=set(pd.read_csv(f"{SP}/T_GCST009379_variants.csv").variant_id)
vx=set(pd.read_csv(f"{SP}/t2d_variants.csv").variant_id)
print(f"\nvariant overlap: Mahajan {len(vm)}, Xue {len(vx)}, shared {len(vm&vx)} ({len(vm&vx)/min(len(vm),len(vx))*100:.0f}% of smaller)")
