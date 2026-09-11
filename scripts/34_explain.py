import sys,numpy as np,pandas as pd
from scipy.stats import spearmanr
SP=sys.argv[1]
MARKS=['H3K27ac','H3K36me3','H3K4me1','H3K4me3','H3K9ac']
h=pd.read_csv(f"{SP}/tracks_chip_histone.csv")
GR=['ISLET','LIVER','ADIPOSE','MUSCLE','BRAIN','IMMUNE','HEART','ARTERY','LUNG','INTESTINE','SKIN']
ok=h['histone_mark'].isin(MARKS).values
IDX={g:((h['grp']==g).values&ok) for g in GR}
NTR={g:int(IDX[g].sum()) for g in GR}
acc=pd.read_csv(f"{SP}/tracks_atac.csv")['grp'].value_counts().add(pd.read_csv(f"{SP}/tracks_dnase.csv")['grp'].value_counts(),fill_value=0).fillna(0)
BG=np.abs(np.load(f"{SP}/bg_scores.npz")['hist'])
R=pd.read_csv(f"{SP}/bench_summary.csv")
# visibility computed from stored matrices
M=pd.read_csv(f"{SP}/bench_matrix.csv",index_col=0)
R2=R.set_index('trait')
vis=M.max(axis=1)                      # strongest tissue signal the trait produces anywhere
R2['visibility']=vis
R2['n_tracks_exp']=[NTR[e] for e in R2['expected']]
R2['acc_tracks_exp']=[int(acc.get(e,0) or 0) for e in R2['expected']]
R2['hit1']=(R2['rank']==1).astype(int)

print("=== Does recovery track ASSAY DEPTH of the expected tissue? ===")
for col in ['n_tracks_exp','acc_tracks_exp']:
    rho,p=spearmanr(R2[col],R2['exp_z']); print(f"  spearman({col}, exp_z) rho={rho:+.3f} p={p:.3f}")
    rho,p=spearmanr(R2[col],-R2['rank']); print(f"  spearman({col}, -rank)  rho={rho:+.3f} p={p:.3f}")

print("\n=== Does recovery track how VISIBLE the trait's variants are to the model? ===")
rho,p=spearmanr(R2['visibility'],R2['exp_z']); print(f"  spearman(visibility, exp_z) rho={rho:+.3f} p={p:.4f}")
rho,p=spearmanr(R2['visibility'],-R2['rank']); print(f"  spearman(visibility, -rank) rho={rho:+.3f} p={p:.4f}")

print("\n=== brain: well instrumented, still fails ===")
print(f"  BRAIN histone tracks={NTR['BRAIN']}, accessibility={int(acc.get('BRAIN',0) or 0)}")
print(R2.loc[[t for t in R2.index if R2.loc[t,'expected']=='BRAIN'],['rank','exp_z','p_exp','visibility']].to_string())
print(f"  ISLET histone tracks={NTR['ISLET']}, accessibility={int(acc.get('ISLET',0) or 0)}")
print(R2.loc[[t for t in R2.index if R2.loc[t,'expected']=='ISLET'],['rank','exp_z','p_exp','visibility']].to_string())

print("\n=== recovered vs not: visibility contrast ===")
a=R2[R2.q_exp<0.05]['visibility']; b=R2[R2.q_exp>=0.05]['visibility']
print(f"  significant traits (n={len(a)}): median visibility {a.median():.2f}")
print(f"  non-significant   (n={len(b)}): median visibility {b.median():.2f}")
R2.to_csv(f"{SP}/bench_explain.csv")
print("\ntrack counts per group:",NTR)
