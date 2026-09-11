import sys,numpy as np,pandas as pd
from scipy.stats import mannwhitneyu,wilcoxon
SP=sys.argv[1]; rng=np.random.default_rng(0)
hist_t=pd.read_csv(f"{SP}/hist_tracks_grp.csv")
MARKS=['H3K27ac','H3K36me3','H3K4me1','H3K4me3']
hok=hist_t['histone_mark'].isin(MARKS).values
GRPS=['ISLET','PANCREAS','LIVER','ADIPOSE','MUSCLE']
IDX={g:((hist_t['grp']==g).values&hok) for g in GRPS}
def raw(t):
    H=np.abs(np.load(f"{SP}/{t}_scores.npz")['hist'])
    return {g:H[:,IDX[g]].mean(1) for g in GRPS}
BG=raw('bg'); NR={g:(BG[g].mean(),BG[g].std()) for g in GRPS}
S={t:pd.DataFrame({g:(raw(t)[g]-NR[g][0])/NR[g][1] for g in GRPS}) for t in ('t2d','ldl','bg')}
def boot(x,n=4000):
    b=rng.choice(x,(n,len(x)),replace=True).mean(1); return np.percentile(b,[2.5,97.5])
print("="*76); print("HISTONE-ONLY (no gene-selection step -> free of that bias)"); print("="*76)
for t in ('t2d','ldl'):
    print(f"\n--- {t.upper()} (n={len(S[t])}) ---")
    r=[]
    for g in GRPS:
        x=S[t][g].values; lo,hi=boot(x)
        r.append(dict(tissue=g,mean_z=round(x.mean(),3),median_z=round(np.median(x),3),
            ci=f"[{lo:.2f},{hi:.2f}]", frac_z_gt1=f"{(x>1).mean()*100:.0f}%",
            p=f"{mannwhitneyu(x,S['bg'][g].values,alternative='greater').pvalue:.1e}"))
    print(pd.DataFrame(r).sort_values('mean_z',ascending=False).to_string(index=False))
print("\n--- within-T2D paired: ISLET vs others ---")
for g in [x for x in GRPS if x!='ISLET']:
    print(f"   ISLET vs {g:9s} diff={S['t2d']['ISLET'].mean()-S['t2d'][g].mean():+.3f}  p={wilcoxon(S['t2d']['ISLET'],S['t2d'][g]).pvalue:.2e}")
print("\n--- within-LDL paired: LIVER vs others ---")
for g in [x for x in GRPS if x!='LIVER']:
    print(f"   LIVER vs {g:9s} diff={S['ldl']['LIVER'].mean()-S['ldl'][g].mean():+.3f}  p={wilcoxon(S['ldl']['LIVER'],S['ldl'][g]).pvalue:.2e}")
for t in ('t2d','ldl'):
    m=pd.read_csv(f"{SP}/{t}_scored_meta.csv").reset_index(drop=True)
    pd.concat([m,S[t].add_suffix('_hz')],axis=1).to_csv(f"{SP}/{t}_histone_z.csv",index=False)
print("\nbackground check:",S['bg'][GRPS].mean().round(3).to_dict())
