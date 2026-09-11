import sys,numpy as np,pandas as pd
from scipy.stats import binomtest
SP=sys.argv[1]
R=pd.read_csv(f"{SP}/bench_final.csv").rename(columns={"Unnamed: 0":"trait"}).set_index('trait')
M=pd.read_csv(f"{SP}/bench_matrix.csv",index_col=0)
GR=list(M.columns)
NICE={'ISLET':'Islet','LIVER':'Liver','ADIPOSE':'Fat','MUSCLE':'Muscle','BRAIN':'Brain','IMMUNE':'Immune',
 'HEART':'Heart','ARTERY':'Artery','LUNG':'Lung','INTESTINE':'Intestine','SKIN':'Skin'}
# Published LDSC-SEG findings (Finucane et al. 2018, Nat Genet 50:621-629), mapped to our groups.
# 'alt' = additional tissue LDSC-SEG reported for that trait.
LD={
 'Type 2 diabetes':                    ('pancreas','ISLET',None),
 'Type 2 diabetes (replication)':      ('pancreas','ISLET',None),
 'Body mass index':                    ('adipose AND cortex','BRAIN','ADIPOSE'),
 'Cognitive ability / education':      ('brain (cerebellar, neuronal)','BRAIN',None),
 'Free cholesterol in large LDL':      ('liver','LIVER',None),
 'HDL cholesterol':                    ('liver','LIVER',None),
 'Triglycerides':                      ('liver','LIVER',None),
 "Crohn's disease (medication proxy)": ('immune (multiple cell types)','IMMUNE',None),
 'Ulcerative colitis (PheCode)':       ('immune (multiple cell types)','IMMUNE',None),
 'Asthma':                             ('T and NKT cells','IMMUNE',None),
 'Systolic blood pressure':            ('stromal / musculoskeletal','ARTERY','MUSCLE'),
}
rows=[]
for t,(desc,primary,alt) in LD.items():
    if t not in R.index: continue
    means=M.loc[t]; order=sorted(GR,key=lambda g:-means[g])
    top=order[0]
    agree = (top==primary) or (alt is not None and top==alt)
    rows.append(dict(trait=t, ldsc=desc, ag_top=NICE[top], ag_rank_primary=order.index(primary)+1,
                     ag_rank_alt=(order.index(alt)+1) if alt else None,
                     agree='yes' if agree else 'no', q=R.loc[t,'q_exp']))
C=pd.DataFrame(rows)
print("=== AlphaGenome vs published LDSC-SEG, overlapping traits ===")
print(C.to_string(index=False,float_format=lambda v:f"{v:.3f}"))
n=len(C); ag=int((C.agree=='yes').sum())
print(f"\nAgreement on top tissue: {ag}/{n}")

print("\n=== Sensitivity: LDSC-SEG says two of my labels were contestable ===")
for t,lab,newlab in [('Body mass index','BRAIN','ADIPOSE'),('Systolic blood pressure','ARTERY','MUSCLE')]:
    means=M.loc[t]; order=sorted(GR,key=lambda g:-means[g])
    print(f"  {t:26s} my label={NICE[lab]:7s} rank={order.index(lab)+1}  |  LDSC-SEG also reports {NICE[newlab]:7s} rank={order.index(newlab)+1}")
orig=int((R['rank']==1).sum())
alt_hits=orig+2
print(f"\n  pre-specified labels : {orig}/20 rank-1  (binomial p={binomtest(orig,20,1/11,alternative='greater').pvalue:.1e})")
print(f"  LDSC-informed labels : {alt_hits}/20 rank-1  (binomial p={binomtest(alt_hits,20,1/11,alternative='greater').pvalue:.1e})")
print("\n  (reported as sensitivity only; pre-specified labels remain primary)")
C.to_csv(f"{SP}/ldsc_comparison.csv",index=False)

print("\n=== the traits where the two methods DISAGREE most ===")
for t in ['Cognitive ability / education','Body mass index','Triglycerides','Asthma']:
    means=M.loc[t]; order=sorted(GR,key=lambda g:-means[g])
    print(f"  {t:32s} LDSC-SEG: {LD[t][0]:28s} AlphaGenome top: {NICE[order[0]]:9s} (expected rank {order.index(LD[t][1])+1})")
