import sys,numpy as np,pandas as pd
from scipy.stats import mannwhitneyu
SP=sys.argv[1]
rna_t=pd.read_csv(f"{SP}/rna_tracks_grp.csv"); hist_t=pd.read_csv(f"{SP}/hist_tracks_grp.csv")
MARKS=['H3K27ac','H3K36me3','H3K4me1','H3K4me3']
hist_ok = hist_t['histone_mark'].isin(MARKS).values
GRPS=['ISLET','PANCREAS','LIVER','ADIPOSE','MUSCLE']

d={}
for tag in ('t2d','ldl'):
    z=np.load(f"{SP}/{tag}_scores.npz")
    d[tag]=dict(rna=np.abs(z['rna']), hist=np.abs(z['hist']),
                meta=pd.read_csv(f"{SP}/{tag}_scored_meta.csv"))
    d[tag]['meta']['trait']=tag

# pooled per-track z-scoring (common reference across both traits)
def zpool(key):
    A=np.vstack([d['t2d'][key], d['ldl'][key]])
    mu,sd=A.mean(0),A.std(0); sd[sd==0]=1
    n=len(d['t2d'][key])
    Z=(A-mu)/sd
    return Z[:n], Z[n:]
rz_t2d,rz_ldl = zpool('rna'); hz_t2d,hz_ldl = zpool('hist')

def tissue_scores(rz,hz):
    out={}
    for g in GRPS:
        ri=(rna_t['grp']==g).values; hi=((hist_t['grp']==g).values & hist_ok)
        parts=[rz[:,ri].mean(1), hz[:,hi].mean(1)]
        out[g]=np.mean(parts,axis=0)          # RNA and histone weighted equally
        out[g+'_rna']=parts[0]; out[g+'_hist']=parts[1]
    return pd.DataFrame(out)
T=tissue_scores(rz_t2d,hz_t2d); L=tissue_scores(rz_ldl,hz_ldl)
T=pd.concat([d['t2d']['meta'].reset_index(drop=True),T],axis=1)
L=pd.concat([d['ldl']['meta'].reset_index(drop=True),L],axis=1)
T.to_csv(f"{SP}/t2d_tissue_scores.csv",index=False); L.to_csv(f"{SP}/ldl_tissue_scores.csv",index=False)

print("="*68); print("MEAN TISSUE SCORE (pooled-z; higher = variant hits this tissue harder)"); print("="*68)
tab=pd.DataFrame({'T2D':[T[g].mean() for g in GRPS],'LDL':[L[g].mean() for g in GRPS]},index=GRPS).round(3)
tab['T2D-LDL']=(tab['T2D']-tab['LDL']).round(3)
tab['MWU_p']=[f"{mannwhitneyu(T[g],L[g],alternative='two-sided').pvalue:.2e}" for g in GRPS]
print(tab.to_string())

print("\n--- CONTROL CHECK: does LDL rank LIVER first? ---")
print(tab['LDL'].sort_values(ascending=False).round(3).to_string())
print("\n--- T2D ranking ---")
print(tab['T2D'].sort_values(ascending=False).round(3).to_string())

print("\n=== per-locus winning tissue ===")
for nm,D in (('T2D',T),('LDL',L)):
    w=D[GRPS].idxmax(axis=1).value_counts()
    print(f"\n{nm} (n={len(D)}):"); print((w/len(D)*100).round(1).astype(str).add('%').to_string())
