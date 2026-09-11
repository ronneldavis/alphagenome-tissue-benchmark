import sys,glob,os,numpy as np,pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
SP=sys.argv[1]
MARKS=['H3K27ac','H3K36me3','H3K4me1','H3K4me3','H3K9ac']
h=pd.read_csv(f"{SP}/tracks_chip_histone.csv"); ok=h['histone_mark'].isin(MARKS).values
GR=['ISLET','LIVER','ADIPOSE','MUSCLE','BRAIN','IMMUNE','HEART','ARTERY','LUNG','INTESTINE','SKIN']
IDX={g:((h['grp']==g).values&ok) for g in GR}
ALL=np.zeros(len(h),bool)
for g in GR: ALL|=IDX[g]
BG=np.abs(np.load(f"{SP}/bg_scores.npz")['hist'])
bgv=BG[:,ALL].mean(1); mu,sd=bgv.mean(),bgv.std()          # tissue-agnostic null

FILES={'Type 2 diabetes':'T_GCST009379','Free cholesterol in large LDL':'T_GCST90501150',
'HDL cholesterol':'T_GCST90501112','Triglycerides':'T_GCST006613','Alanine aminotransferase':'T_GCST90474304',
'Alkaline phosphatase':'T_GCST90018942',"Crohn's disease (medication proxy)":'T_FINNGEN_R12_RX_CROHN_1STLINE',
'Ulcerative colitis (PheCode)':'T_GCST90480318','Type 1 diabetes (PheCode)':'T_GCST90479877',
'Asthma':'T_GCST010043','Platelet count':'T_GCST90662907','Atrial fibrillation':'T_GCST90624411',
'Heart rate':'T_GCST90480666','Coronary artery disease':'T_GCST90132314','Systolic blood pressure':'T_GCST90000066',
'Body mass index':'T_GCST007039','Cognitive ability / education':'T_GCST008595',
'Lung function (FEV1/FVC)':'T_GCST90244094','Atopic dermatitis':'T_GCST90503108',
'Type 2 diabetes (replication)':'t2d'}
R=pd.read_csv(f"{SP}/bench_summary.csv").set_index('trait')
vis={}
for t,f in FILES.items():
    H=np.abs(np.load(f"{SP}/{f}_scores.npz")['hist'])
    vis[t]=((H[:,ALL].mean(1)-mu)/sd).mean()
R['visibility']=pd.Series(vis)
R['recovered']=(R['rank']==1).astype(int)
E=pd.read_csv(f"{SP}/bench_explain.csv",index_col=0)
R['n_tracks_exp']=E['n_tracks_exp']; R['acc_tracks_exp']=E['acc_tracks_exp']

print("=== TISSUE-AGNOSTIC visibility (independent of which tissue is expected) ===")
for a,b,lab in [('visibility','exp_z','expected-tissue enrichment'),
                ('n_tracks_exp','exp_z','expected-tissue enrichment'),
                ('acc_tracks_exp','exp_z','expected-tissue enrichment')]:
    rho,p=spearmanr(R[a],R[b]); print(f"  spearman({a:15s}, {lab}) rho={rho:+.3f}  p={p:.4f}")
rho,p=spearmanr(R['visibility'],-R['rank']); print(f"  spearman(visibility     , -rank) rho={rho:+.3f}  p={p:.4f}")
a=R[R.q_exp<0.05]['visibility']; b=R[R.q_exp>=0.05]['visibility']
u=mannwhitneyu(a,b,alternative='greater')
print(f"\n  recovered (q<0.05, n={len(a)}) median visibility {a.median():+.3f}")
print(f"  not       (q>=0.05, n={len(b)}) median visibility {b.median():+.3f}   MWU p={u.pvalue:.4f}")
print("\n=== the decisive contrast ===")
for t in ['Type 2 diabetes','Type 2 diabetes (replication)','Body mass index','Cognitive ability / education',
          'Alanine aminotransferase','Platelet count']:
    r=R.loc[t]; print(f"  {t:32s} exp={r['expected']:7s} tracks={int(r['n_tracks_exp']):3d} acc={int(r['acc_tracks_exp']):3d} "
                      f"rank={int(r['rank'])} exp_z={r['exp_z']:+.2f} vis={r['visibility']:+.2f}")
R.to_csv(f"{SP}/bench_final.csv")
