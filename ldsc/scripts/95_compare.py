import os,sys,json,numpy as np,pandas as pd
from scipy.stats import binomtest
L=sys.argv[1]; SP=os.path.dirname(L)
U=json.load(open(f"{L}/sumstats_urls.json"))
AG=pd.read_csv(f"{SP}/bench_final.csv").rename(columns={"Unnamed: 0":"trait"}).set_index("trait")
AGM=pd.read_csv(f"{SP}/bench_matrix.csv",index_col=0)
GR=list(AGM.columns)
NICE={'ISLET':'Islet','LIVER':'Liver','ADIPOSE':'Fat','MUSCLE':'Muscle','BRAIN':'Brain','IMMUNE':'Immune',
 'HEART':'Heart','ARTERY':'Artery','LUNG':'Lung','INTESTINE':'Intestine','SKIN':'Skin'}
rows=[]
for a,d in U.items():
    p=f"{L}/results/{a}.cell_type_results.txt"
    if not os.path.exists(p): continue
    t=pd.read_csv(p,sep='\t').sort_values('Coefficient_P_value').reset_index(drop=True)
    exp=d['tissue']; tr=d['trait']
    ld_rank=list(t.Name).index(exp)+1
    ld_top=t.Name[0]; ld_p=float(t.loc[t.Name==exp,'Coefficient_P_value'].iloc[0])
    if tr in AGM.index:
        m=AGM.loc[tr]; order=sorted(GR,key=lambda g:-m[g])
        ag_rank=order.index(exp)+1; ag_top=order[0]
        ag_q=float(AG.loc[tr,'q_exp'])
    else:
        ag_rank=np.nan; ag_top='--'; ag_q=np.nan
    rows.append(dict(trait=tr,expected=exp,ag_rank=ag_rank,ag_top=ag_top,ag_q=ag_q,
                     ld_rank=ld_rank,ld_top=ld_top,ld_p=ld_p,
                     agree=(ag_top==ld_top)))
C=pd.DataFrame(rows)
C.to_csv(f"{SP}/method_comparison.csv",index=False)
n=len(C)
print("="*104)
print("MATCHED HEAD-TO-HEAD: AlphaGenome vs LDSC-SEG  (same traits, same GWAS, same 11 tissue definitions)")
print("="*104)
pr=C.copy()
pr['expected']=pr.expected.map(NICE); pr['ag_top']=pr.ag_top.map(lambda x:NICE.get(x,x)); pr['ld_top']=pr.ld_top.map(NICE)
pr=pr[['trait','expected','ag_rank','ag_top','ld_rank','ld_top','ld_p','agree']]
pr.columns=['Trait','Expected','AG rank','AG top','LDSC rank','LDSC top','LDSC P(exp)','Agree']
print(pr.sort_values('LDSC rank').to_string(index=False,float_format=lambda v:f"{v:.3g}"))
ag1=int((C.ag_rank==1).sum()); ld1=int((C.ld_rank==1).sum())
ag3=int((C.ag_rank<=3).sum()); ld3=int((C.ld_rank<=3).sum())
print(f"\n  rank-1 recovery : AlphaGenome {ag1}/{n}   LDSC-SEG {ld1}/{n}   (chance {n/11:.1f})")
print(f"  top-3 recovery  : AlphaGenome {ag3}/{n}   LDSC-SEG {ld3}/{n}")
print(f"  binomial vs 1/11: AG P={binomtest(ag1,n,1/11,alternative='greater').pvalue:.2e}  "
      f"LDSC P={binomtest(ld1,n,1/11,alternative='greater').pvalue:.2e}")
print(f"  methods agree on top tissue: {int(C.agree.sum())}/{n}")
both=C[(C.ag_rank>1)&(C.ld_rank>1)]
print(f"\n  traits BOTH methods miss: {len(both)}  -> {list(both.trait)}")
print(f"  of those, both pick the SAME wrong tissue: {int((both.ag_top==both.ld_top).sum())}")
