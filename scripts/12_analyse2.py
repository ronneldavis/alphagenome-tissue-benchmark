import sys,numpy as np,pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon
SP=sys.argv[1]; rng=np.random.default_rng(0)
rna_t=pd.read_csv(f"{SP}/rna_tracks_grp.csv"); hist_t=pd.read_csv(f"{SP}/hist_tracks_grp.csv")
MARKS=['H3K27ac','H3K36me3','H3K4me1','H3K4me3']
hok=hist_t['histone_mark'].isin(MARKS).values
GRPS=['ISLET','PANCREAS','LIVER','ADIPOSE','MUSCLE']
IDX={g:((rna_t['grp']==g).values,((hist_t['grp']==g).values&hok)) for g in GRPS}

def raw(tag):
    z=np.load(f"{SP}/{tag}_scores.npz"); R,H=np.abs(z['rna']),np.abs(z['hist'])
    return {g:(R[:,IDX[g][0]].mean(1), H[:,IDX[g][1]].mean(1)) for g in GRPS}
BG=raw('bg'); norm={g:((BG[g][0].mean(),BG[g][0].std()),(BG[g][1].mean(),BG[g][1].std())) for g in GRPS}

def score(tag):
    R=raw(tag); out={}
    for g in GRPS:
        (rm,rs),(hm,hs)=norm[g]
        out[g+'_rna']=(R[g][0]-rm)/rs; out[g+'_hist']=(R[g][1]-hm)/hs
        out[g]=(out[g+'_rna']+out[g+'_hist'])/2
    return pd.DataFrame(out)
S={t:score(t) for t in ('t2d','ldl','bg')}
for t in ('t2d','ldl'):
    m=pd.read_csv(f"{SP}/{t}_scored_meta.csv").reset_index(drop=True)
    pd.concat([m,S[t]],axis=1).to_csv(f"{SP}/{t}_tissue_z.csv",index=False)

def boot(x,n=4000):
    b=rng.choice(x,(n,len(x)),replace=True).mean(1); return np.percentile(b,[2.5,97.5])
print("="*84)
print("ENRICHMENT vs MATCHED BACKGROUND  (z>0 = trait variants hit this tissue harder than random)")
print("="*84)
rows=[]
for t in ('t2d','ldl'):
    for g in GRPS:
        x=S[t][g].values; lo,hi=boot(x)
        p=mannwhitneyu(x,S['bg'][g].values,alternative='greater').pvalue
        rows.append(dict(trait=t.upper(),tissue=g,mean_z=round(x.mean(),3),
                         ci=f"[{lo:.2f},{hi:.2f}]",p_vs_bg=f"{p:.1e}",n=len(x)))
R=pd.DataFrame(rows)
for t in ('T2D','LDL'):
    print(f"\n--- {t} ---")
    print(R[R.trait==t].drop(columns='trait').sort_values('mean_z',ascending=False).to_string(index=False))

print("\n"+"="*84); print("WITHIN-TRAIT PAIRED TEST: is the top tissue really above the others?"); print("="*84)
for t in ('t2d','ldl'):
    order=S[t][GRPS].mean().sort_values(ascending=False)
    top=order.index[0]
    print(f"\n{t.upper()}: top tissue = {top}")
    for g in order.index[1:]:
        st=wilcoxon(S[t][top],S[t][g]); d=S[t][top].mean()-S[t][g].mean()
        print(f"   {top} vs {g:9s} diff={d:+.3f}  wilcoxon p={st.pvalue:.2e}")
print("\n"+"="*84); print("SANITY: background should sit at ~0"); print("="*84)
print(S['bg'][GRPS].mean().round(3).to_string())
