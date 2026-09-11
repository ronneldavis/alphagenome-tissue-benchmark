import sys,numpy as np,pandas as pd
from scipy.stats import binomtest, spearmanr, mannwhitneyu
SP=sys.argv[1]; P=f"{SP}/paper"
R=pd.read_csv(f"{SP}/bench_final.csv").rename(columns={"Unnamed: 0":"trait"})
tc={m:pd.read_csv(f"{SP}/tracks_{m}.csv")['grp'].value_counts() for m in ['rna_seq','atac','dnase','chip_histone']}
TC=pd.DataFrame(tc).fillna(0).astype(int)
nbg=int(np.load(f"{SP}/bg_scores.npz")['hist'].shape[0])
n=len(R); NT=11
h1=int((R['rank']==1).sum()); h3=int((R['rank']<=3).sum()); nsig=int((R.q_exp<0.05).sum())
def sci(p):
    if p>=0.01: return f"{p:.2f}"
    m,e=f"{p:.1e}".split("e"); return f"{m}\\times 10^{{{int(e)}}}"
rho_v,p_v=spearmanr(R['visibility'],R['exp_z'])
rho_t,p_t=spearmanr(R['n_tracks_exp'],R['exp_z'])
rho_a,p_a=spearmanr(R['acc_tracks_exp'],R['exp_z'])
a=R[R.q_exp<0.05]['visibility']; b=R[R.q_exp>=0.05]['visibility']
pmw=mannwhitneyu(a,b,alternative='greater').pvalue
g=R.set_index('trait')
M={'NTRAITS':n,'NBG':nbg,'NTISSUE':NT,'NHITONE':h1,'NHITTHREE':h3,'NSIG':nsig,'NCAP':150,
 'HITONEPCT':f"{h1/n*100:.0f}",'HITTHREEPCT':f"{h3/n*100:.0f}",
 'PHITONE':sci(binomtest(h1,n,1/NT,alternative='greater').pvalue),
 'PHITTHREE':sci(binomtest(h3,n,3/NT,alternative='greater').pvalue),
 'CHANCEONE':f"{n/NT:.1f}",'CHANCETHREE':f"{3*n/NT:.1f}",
 'RHOVIS':f"{rho_v:+.2f}",'PVIS':sci(p_v),'RHOTRACK':f"{rho_t:+.2f}",'PTRACK':sci(p_t),
 'RHOACC':f"{rho_a:+.2f}",'PACC':sci(p_a),
 'VISREC':f"{a.median():.2f}",'VISNOT':f"{b.median():.2f}",'PVISMWU':sci(pmw),
 'BRAINHIST':int(TC.loc['BRAIN','chip_histone']),
 'BRAINACC':int(TC.loc['BRAIN','atac']+TC.loc['BRAIN','dnase']),
 'ISLETRNA':int(TC.loc['ISLET','rna_seq']),'ISLETHIST':int(TC.loc['ISLET','chip_histone']),
 'ISLETACC':int(TC.loc['ISLET','atac']+TC.loc['ISLET','dnase']),
 'LIVERACC':int(TC.loc['LIVER','atac']+TC.loc['LIVER','dnase']),
 'IMMUNEACC':int(TC.loc['IMMUNE','atac']+TC.loc['IMMUNE','dnase']),
 'TDZ':f"{g.loc['Type 2 diabetes','exp_z']:.2f}",'TDP':f"{g.loc['Type 2 diabetes','p_exp']:.3f}",
 'TDXZ':f"{g.loc['Type 2 diabetes (replication)','exp_z']:.2f}",
 'TDXRANK':int(g.loc['Type 2 diabetes (replication)','rank']),
 'TDXP':f"{g.loc['Type 2 diabetes (replication)','p_exp']:.2f}",
 'BMIZ':f"{g.loc['Body mass index','exp_z']:.2f}",'BMIRANK':int(g.loc['Body mass index','rank']),
 'COGZ':f"{g.loc['Cognitive ability / education','exp_z']:.2f}",
 'ALTZ':f"{g.loc['Alanine aminotransferase','exp_z']:.2f}",
 'CODEURL':"the repository listed under Data availability",
}
vm=set(pd.read_csv(f"{SP}/T_GCST009379_variants.csv").variant_id)
vx=set(pd.read_csv(f"{SP}/t2d_variants.csv").variant_id)
M['OVERLAP']=len(vm&vx); M['OVERLAPPCT']=f"{len(vm&vx)/min(len(vm),len(vx))*100:.0f}"
M['CSMEDIAN']=int(pd.read_csv(f"{SP}/t2d_variants.csv")['cs_size'].median())
open(f"{P}/numbers.tex","w").write("\n".join(f"\\newcommand{{\\{k}}}{{{v}}}" for k,v in M.items())+"\n")
for k in ['NHITONE','NHITTHREE','NSIG','PHITONE','RHOVIS','PVIS','RHOTRACK','PTRACK','VISREC','VISNOT','PVISMWU','OVERLAP','TDZ','TDP','TDXRANK']:
    print(f"  {k} = {M[k]}")
