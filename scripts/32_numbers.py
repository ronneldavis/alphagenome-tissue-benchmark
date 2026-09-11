import sys,json,numpy as np,pandas as pd
from scipy.stats import binomtest
SP=sys.argv[1]; P=f"{SP}/paper"
R=pd.read_csv(f"{SP}/bench_final.csv").rename(columns={"Unnamed: 0":"trait"})
tc={m:pd.read_csv(f"{SP}/tracks_{m}.csv")['grp'].value_counts() for m in ['rna_seq','atac','dnase','chip_histone']}
TC=pd.DataFrame(tc).fillna(0).astype(int)
nbg=int(np.load(f"{SP}/bg_scores.npz")['hist'].shape[0])
n=len(R); NT=11
h1=int((R['rank']==1).sum()); h3=int((R['rank']<=3).sum())
p1=binomtest(h1,n,1/NT,alternative='greater').pvalue
p3=binomtest(h3,n,3/NT,alternative='greater').pvalue
t2=R[R.trait.str.startswith('Type 2')]
ranks=sorted(t2['rank'].tolist())
rk=" and ".join(str(x) for x in ranks) if len(ranks)>1 else str(ranks[0] if ranks else "n/a")
def fm(p): return f"{p:.1e}".replace("e-0","\\times 10^{-").replace("e-","\\times 10^{-")+"}" if p<0.01 else f"{p:.2f}"
cs=[]
for f in ["t2d_variants.csv"]:
    try: cs.append(pd.read_csv(f"{SP}/{f}")['cs_size'].median())
    except Exception: pass
M={
 'NTRAITS':n,'NBG':nbg,'NTISSUE':NT,'NHITONE':h1,'NHITTHREE':h3,
 'HITONEPCT':f"{h1/n*100:.0f}",'HITTHREEPCT':f"{h3/n*100:.0f}",
 'PHITONE':fm(p1),'PHITTHREE':fm(p3),'CHANCEONE':f"{n/NT:.1f}",
 'TWODRANK':rk,'NCAP':150,
 'ISLETRNA':TC.loc['ISLET','rna_seq'],'ISLETHIST':TC.loc['ISLET','chip_histone'],
 'ISLETDNASE':TC.loc['ISLET','dnase'],
 'LIVERACC':int(TC.loc['LIVER','atac']+TC.loc['LIVER','dnase']),
 'IMMUNEACC':int(TC.loc['IMMUNE','atac']+TC.loc['IMMUNE','dnase']),
 'CSMEDIAN':f"the median 95\\% credible set in the primary T2D dataset contains {int(cs[0]) if cs else 5} variants.",
 'CODEURL':"\\url{https://github.com/} (see Data availability)",
}
open(f"{P}/numbers.tex","w").write("\n".join(f"\\newcommand{{\\{k}}}{{{v}}}" for k,v in M.items())+"\n")
print(json.dumps({k:str(v) for k,v in M.items() if k not in ('CSMEDIAN','CODEURL')},indent=0))
