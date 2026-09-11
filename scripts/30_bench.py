import os,sys,json,numpy as np,pandas as pd
from scipy.stats import binomtest, wilcoxon, mannwhitneyu
SP=sys.argv[1]
MARKS=['H3K27ac','H3K36me3','H3K4me1','H3K4me3','H3K9ac']
h=pd.read_csv(f"{SP}/tracks_chip_histone.csv")
GR=[g for g in ['ISLET','LIVER','ADIPOSE','MUSCLE','BRAIN','IMMUNE','HEART','ARTERY','LUNG','INTESTINE','SKIN']]
ok=h['histone_mark'].isin(MARKS).values
IDX={g:((h['grp']==g).values & ok) for g in GR}
BG=np.abs(np.load(f"{SP}/bg_scores.npz")['hist'])
NR={g:(BG[:,IDX[g]].mean(1).mean(), BG[:,IDX[g]].mean(1).std()) for g in GR}

LABEL={ # study -> (display trait, expected tissue, note)
 'GCST009379':('Type 2 diabetes','ISLET','Mahajan 2018'),
 'GCST006867':('Type 2 diabetes (replication)','ISLET','Xue 2018'),
 'GCST90501150':('Free cholesterol in large LDL','LIVER','Zoodsma 2025'),
 'GCST90501112':('HDL cholesterol','LIVER','Zoodsma 2025'),
 'GCST006613':('Triglycerides','LIVER','Klarin 2018'),
 'GCST90474304':('Alanine aminotransferase','LIVER','UKB-WGS 2025'),
 'GCST90018942':('Alkaline phosphatase','LIVER','Sakaue 2021'),
 'FINNGEN_R12_RX_CROHN_1STLINE':("Crohn's disease (medication proxy)",'IMMUNE','FinnGen R12'),
 'GCST90480318':('Ulcerative colitis (PheCode)','IMMUNE','Verma 2024'),
 'GCST90479877':('Type 1 diabetes (PheCode)','IMMUNE','Verma 2024'),
 'GCST010043':('Asthma','IMMUNE','Han 2020'),
 'GCST90662907':('Platelet count','IMMUNE','Jee 2025'),
 'GCST90624411':('Atrial fibrillation','HEART','Yuan 2025'),
 'GCST90480666':('Heart rate','HEART','Verma 2024'),
 'GCST90132314':('Coronary artery disease','ARTERY','Aragam 2022'),
 'GCST90000066':('Systolic blood pressure','ARTERY','Surendran 2020'),
 'GCST007039':('Body mass index','BRAIN','Kichaev 2018'),
 'GCST008595':('Cognitive ability / education','BRAIN','Lam 2019'),
 'GCST90244094':('Lung function (FEV1/FVC)','LUNG','Shrine 2023'),
 'GCST90503108':('Atopic dermatitis','SKIN','Oliva 2025'),
}
rows=[];mat={};nvar={}
for sid,(name,exp,src) in LABEL.items():
    f=f"{SP}/T_{sid}_scores.npz"
    if not os.path.exists(f): continue
    H=np.abs(np.load(f)['hist'])
    if H.ndim!=2 or H.shape[0]<20: continue
    z={g:((H[:,IDX[g]].mean(1)-NR[g][0])/NR[g][1]) for g in GR}
    means={g:z[g].mean() for g in GR}
    order=sorted(GR,key=lambda g:-means[g])
    rank=order.index(exp)+1
    mat[name]=means; nvar[name]=H.shape[0]
    rows.append(dict(trait=name,source=src,n=H.shape[0],expected=exp,
      rank=rank,top=order[0],top_z=round(means[order[0]],3),exp_z=round(means[exp],3),
      p_exp=mannwhitneyu(z[exp],(BG[:,IDX[exp]].mean(1)-NR[exp][0])/NR[exp][1],alternative='greater').pvalue))
R=pd.DataFrame(rows).sort_values('rank')
print(R.to_string(index=False,float_format=lambda v:f"{v:.2e}" if v<0.01 else f"{v:.3f}"))
n=len(R); hit1=int((R['rank']==1).sum()); hit3=int((R['rank']<=3).sum())
print(f"\nTraits: {n} | expected tissue ranked #1: {hit1} ({hit1/n*100:.0f}%) | top-3: {hit3} ({hit3/n*100:.0f}%)")
print(f"chance: #1 = {n/11:.1f} traits, top-3 = {3*n/11:.1f}")
print("binomial p (rank1 vs 1/11):", f"{binomtest(hit1,n,1/11,alternative='greater').pvalue:.2e}")
print("binomial p (top3 vs 3/11):", f"{binomtest(hit3,n,3/11,alternative='greater').pvalue:.2e}")
print("\nmedian rank of expected tissue:",R['rank'].median(),"| T2D ranks:",R[R.trait.str.startswith('Type 2')][['trait','rank','top','exp_z']].to_string(index=False))
pd.DataFrame(mat).T[GR].to_csv(f"{SP}/bench_matrix.csv")
R.to_csv(f"{SP}/bench_summary.csv",index=False)
print("\nsaved bench_matrix.csv / bench_summary.csv")
