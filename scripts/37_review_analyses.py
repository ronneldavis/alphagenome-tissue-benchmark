"""Review-stage analyses added after internal review (Sept 2026).

Usage: python scripts/37_review_analyses.py <SCORE_DIR> <PAPER_DIR> <PANEL_JSON>
  SCORE_DIR holds bg_scores.npz, T_<study>_scores.npz, t2d_scores.npz and the *_variants.csv files
  produced by 22_run_all.py / 11_background.py; PAPER_DIR is this repository. PANEL_JSON is the
  panel definition file (e.g. data/panel_v1.json): a list of {trait, label, study, expected,
  source, file} dicts, used here for the trait->file mapping (FILES) and the trait->display-label
  mapping (LABEL, taken only where label != trait).
Writes: numbers_extra.tex, results_table.tex (with bootstrap column), traits_table.tex (with GWAS
size, case counts and fine-mapping regime), tracks_table.tex (with five-mark histone column),
fig_heatmap.pdf/png (relabelled), ldsc_table.tex (relabelled), data/pervariant_summary.csv.
"""
import sys, os, json, numpy as np, pandas as pd
from scipy.stats import spearmanr, mannwhitneyu, binomtest
import statsmodels.formula.api as smf
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
SP, P, PANEL = sys.argv[1], sys.argv[2], sys.argv[3]; D = f"{P}/data"
PANEL_ROWS = json.load(open(PANEL))
MARKS=['H3K27ac','H3K36me3','H3K4me1','H3K4me3','H3K9ac']
GR=['ISLET','LIVER','ADIPOSE','MUSCLE','BRAIN','IMMUNE','HEART','ARTERY','LUNG','INTESTINE','SKIN']
NICE={'ISLET':'Islet','LIVER':'Liver','ADIPOSE':'Fat','MUSCLE':'Muscle','BRAIN':'Brain','IMMUNE':'Immune',
      'HEART':'Heart','ARTERY':'Artery','LUNG':'Lung','INTESTINE':'Intestine','SKIN':'Skin'}
LABEL={e['trait']:e['label'] for e in PANEL_ROWS if e['label']!=e['trait']}
lab=lambda t: LABEL.get(t,t)
esc=lambda s:str(s).replace('&','\\&').replace('_','\\_').replace('%','\\%')
def sci(p):
    if p>=0.01: return f"{p:.2f}"
    m,e=f"{p:.1e}".split("e"); return f"{m}\\times 10^{{{int(e)}}}"
h=pd.read_csv(f"{D}/tracks_chip_histone.csv"); ok=h['histone_mark'].isin(MARKS).values
IDX={g:((h['grp']==g).values&ok) for g in GR}; ALL=np.zeros(len(h),bool)
for g in GR: ALL|=IDX[g]
FIVE={g:int(IDX[g].sum()) for g in GR}
BG=np.abs(np.load(f"{SP}/bg_scores.npz")['hist'])
NR={g:(BG[:,IDX[g]].mean(1).mean(), BG[:,IDX[g]].mean(1).std()) for g in GR}
ZB=np.column_stack([(BG[:,IDX[g]].mean(1)-NR[g][0])/NR[g][1] for g in GR])
null_rate={g:float((ZB.argmax(1)==i).mean()) for i,g in enumerate(GR)}
bgv=BG[:,ALL].mean(1); mu,sd=bgv.mean(),bgv.std()
FILES={e['trait']:e['file'] for e in PANEL_ROWS}
# NOTE: the bootstrap loop below draws from one shared `rng`, in FILES's
# iteration order (= the panel file's row order), so the *order of rows in
# the panel JSON* is part of this script's reproducibility contract -- do
# not reorder an existing panel file's rows, or every bootstrap proportion
# from that row onward will change (a real fixed-seed effect, not a bug).
# See scripts/README.md.
R=pd.read_csv(f"{D}/bench_final.csv").set_index('trait')
rng=np.random.default_rng(1); B=2000
rows=[]; allv=[]; Zall=[]
for t,f in FILES.items():
    H=np.abs(np.load(f"{SP}/{f}_scores.npz")['hist']); exp=R.loc[t,'expected']; ei=GR.index(exp)
    Z=np.column_stack([(H[:,IDX[g]].mean(1)-NR[g][0])/NR[g][1] for g in GR]); vis=(H[:,ALL].mean(1)-mu)/sd
    vrank=(Z>Z[:,[ei]]).sum(1)+1
    n=len(Z); b1=sum(Z[rng.integers(0,n,n)].mean(0).argmax()==ei for _ in range(B))/B
    var=pd.read_csv(f"{SP}/{f}_variants.csv"); pcol='negp' if 'negp' in var else 'neglog10p'
    if os.path.exists(f"{SP}/{f}_meta.csv"):
        meta=pd.read_csv(f"{SP}/{f}_meta.csv"); mm=meta.merge(var[['variant_id',pcol]],on='variant_id',how='left')
        negp=mm[pcol].values; pip=meta.pip.values
    else:
        var=var[~var.coding].reset_index(drop=True)
        negp=var[pcol].values if len(var)==n else np.full(n,np.nan); pip=var.pip.values if len(var)==n else np.full(n,np.nan)
    rows.append(dict(trait=t,n=n,boot=b1,vis=vis.mean(),med_negp=np.nanmedian(negp)))
    allv.append(pd.DataFrame(dict(trait=t,exp=exp,vis=vis,vrank=vrank,null=null_rate[exp],negp=negp,pip=pip))); Zall.append(Z)
T=pd.DataFrame(rows).set_index('trait'); R=R.join(T[['boot','med_negp']])
V=pd.concat(allv); V['hit']=(V.vrank==1).astype(int)
# ---- trait-level statistics
M={}
M['NUNIQUE']=19; M['NVAR']=f"{len(V):,}"; M['NNOTSIG']=int((R.q_exp>=0.05).sum())
h1=int((R['rank']==1).sum()); M['NFRAGILE']=int(((R['rank']==1)&(R.boot<0.5)).sum())
M['PHITONEUNIQ']=sci(binomtest(h1,19,1/11,alternative='greater').pvalue)
M['PHITONEUNIQB']=sci(binomtest(h1-1,19,1/11,alternative='greater').pvalue)
rho,p=spearmanr(R.visibility,-R['rank']); M['RHOVISRANK']=f"{rho:+.2f}"; M['PVISRANK']=sci(p)
rho,p=spearmanr(R.visibility,R.boot); M['RHOVISBOOT']=f"{rho:+.2f}"; M['PVISBOOT']=sci(p)
rho,p=spearmanr(R.visibility,R.top_z); M['RHOVISTOP']=f"{rho:+.2f}"; M['PVISTOP']=sci(p)
a=R[R['rank']==1].visibility; b=R[R['rank']!=1].visibility
M['VISRANKONE']=f"{a.median():.2f}"; M['VISRANKNOT']=f"{b.median():.2f}"; M['PVISRANKMWU']=sci(mannwhitneyu(a,b,alternative='greater').pvalue)
rho,p=spearmanr(R.n_tracks_exp,-R['rank']); M['RHOTRACKRANK']=f"{rho:+.2f}"; M['PTRACKRANK']=sci(p)
rho,p=spearmanr(R.acc_tracks_exp,-R['rank']); M['RHOACCRANK']=f"{rho:+.2f}"; M['PACCRANK']=sci(p)
rho,_=spearmanr(R.n_tracks_exp,R.exp_z); z=np.arctanh(rho); se=1/np.sqrt(len(R)-3)
M['RHOTRACKLO']=f"{np.tanh(z-1.96*se):+.2f}"; M['RHOTRACKHI']=f"{np.tanh(z+1.96*se):+.2f}"
o=R.sort_values('visibility',ascending=False); M['NTOPVISHIT']=int((o['rank'].iloc[:10]==1).sum()); M['NBOTVISHIT']=int((o['rank'].iloc[-10:]==1).sum())
rho,p=spearmanr(R.med_negp,R.exp_z); M['RHONEGPZ']=f"{rho:+.2f}"; M['PNEGPZ']=sci(p)
rho,p=spearmanr(R.med_negp,R.visibility); M['RHONEGPVIS']=f"{rho:+.2f}"; M['PNEGPVIS']=sci(p)
rho,p=spearmanr(V.negp,V.vis,nan_policy='omit'); M['RHOVNEGP']=f"{rho:+.2f}"; M['PVNEGP']=sci(p)
rho,p=spearmanr(V.pip,V.vis,nan_policy='omit'); M['RHOVPIP']=f"{rho:+.2f}"; M['PVPIP']=sci(p)
M['TDBOOT']=f"{R.loc['Type 2 diabetes','boot']:.2f}"
# ---- shared component
Zp=np.vstack(Zall); CM=pd.DataFrame(Zp,columns=GR).corr(method="spearman").values; off=CM[np.triu_indices(11,1)]
M['INTERTISSUE']=f"{off.mean():.2f}"; M['INTERTISSUEMIN']=f"{off.min():.2f}"; M['INTERTISSUEMAX']=f"{off.max():.2f}"
s=np.linalg.svd(Zp-Zp.mean(0),compute_uv=False); M['PCONE']=f"{100*s[0]**2/(s**2).sum():.0f}"
# ---- per-variant calibrated analysis
M['IMMUNENULL']=f"{100*null_rate['IMMUNE']:.1f}"; M['LIVERNULL']=f"{100*null_rate['LIVER']:.1f}"
M['ALLOBS']=f"{100*V.hit.mean():.1f}"; M['ALLNULL']=f"{100*V.null.mean():.1f}"; M['ALLRATIO']=f"{V.hit.mean()/V.null.mean():.2f}"
V['q']=pd.qcut(V.vis,4,labels=['Q1','Q2','Q3','Q4']); qs={}
for q,sq in V.groupby('q',observed=True):
    qs[q]=(sq.hit.mean()/sq.null.mean(), binomtest(int(sq.hit.sum()),len(sq),sq.null.mean(),alternative='greater').pvalue, 100*sq.hit.mean(), 100*sq.null.mean())
M['QLOWRATIO']=f"{qs['Q1'][0]:.2f}"; M['QTWORATIO']=f"{qs['Q2'][0]:.2f}"; M['QTHREERATIO']=f"{qs['Q3'][0]:.2f}"
M['QHIGHRATIO']=f"{qs['Q4'][0]:.2f}"; M['QHIGHP']=sci(qs['Q4'][1]); M['QHIGHOBS']=f"{qs['Q4'][2]:.1f}"; M['QHIGHNULL']=f"{qs['Q4'][3]:.1f}"
V['visz']=(V.vis-V.vis.mean())/V.vis.std(); V['top']=(V.groupby('trait').vis.rank(pct=True)>0.75).astype(int)
m=smf.logit('hit ~ visz + C(trait)',data=V).fit(disp=0); M['ORVISSD']=f"{np.exp(m.params['visz']):.2f}"; M['PVISSD']=sci(m.pvalues['visz'])
m2=smf.logit('hit ~ top + C(trait)',data=V).fit(disp=0); ci=np.exp(m2.conf_int().loc['top'])
M['ORTOPQ']=f"{np.exp(m2.params['top']):.2f}"; M['ORTOPQLO']=f"{ci[0]:.2f}"; M['ORTOPQHI']=f"{ci[1]:.2f}"; M['PTOPQ']=sci(m2.pvalues['top'])
# ---- calibrated trait-level chance (bootstrap background sets of n=150)
rng2=np.random.default_rng(0); cnt=np.zeros(11)
for _ in range(4000): cnt[ZB[rng2.integers(0,len(ZB),150)].mean(0).argmax()]+=1
tl=dict(zip(GR,cnt/cnt.sum())); mix=R.expected.value_counts().to_dict()
M['CHANCECAL']=f"{sum(v*tl[g] for g,v in mix.items()):.1f}"
# ---- LDSC-derived numbers
C2=pd.read_csv(f"{D}/method_comparison.csv")
M['AGONESIG']=int(((C2.ag_rank==1)&(C2.ag_q<0.05)).sum()); M['LDONESIG']=int(((C2.ld_rank==1)&(C2.ld_p<0.05)).sum())
ld=lambda a: pd.read_csv(f"{P}/ldsc/results/{a}.cell_type_results.txt",sep='\t').set_index('Name')
ad=ld('GCST90503108'); M['ADLDP']=sci(float(ad.loc['IMMUNE','Coefficient_P_value'])); M['ADAGZ']=f"{R.loc['Atopic dermatitis','top_z']:.2f}"
cad=ld('GCST90132314').sort_values('Coefficient_P_value'); M['CADLDTOPP']=sci(float(cad.Coefficient_P_value.iloc[0]))
M['CADLDRANK']=int(list(cad.index).index('ARTERY')+1); M['CADAGQ']=f"{R.loc['Coronary artery disease','q_exp']:.3f}"
# ---- provenance
FM=json.load(open(f"{D}/finemapping_provenance.json")); SM=json.load(open(f"{D}/study_meta.json")); SID={e["trait"]:e["study"] for e in json.load(open(sys.argv[3]))}  # accession per trait, from the panel file
def regime(s):
    conf=FM[s]['confidence']; tot=sum(conf.values())
    top=max(conf,key=conf.get); share=conf[top]/tot
    name={'PICS fine-mapped credible set based on reported top hit':'PICS (top hits)','PICS fine-mapped credible set extracted from summary statistics':'PICS (sumstats)',
          'SuSiE fine-mapped credible set with out-of-sample LD':'SuSiE-inf','SuSiE fine-mapped credible set with in-sample LD':'SuSiE'}[top]
    return name if share>0.99 else f"{name} ({100*share:.0f}\\%)"
reg={t:regime(SID[t]) for t in R.index}
M['NSUSIE']=sum(v.startswith('SuSiE') for v in reg.values()); M['NPICSSUM']=sum(v.startswith('PICS (sumstats') for v in reg.values()); M['NPICSTOP']=sum(v.startswith('PICS (top') for v in reg.values())
M['UCCASES']=f"{FM['GCST90480318']['nCases']:,}"; M['TODCASES']=f"{FM['GCST90479877']['nCases']:,}"; M['CROHNCASES']=f"{FM['FINNGEN_R12_RX_CROHN_1STLINE']['nCases']:,}"
M['BRAINFIVE']=FIVE['BRAIN']; M['ISLETFIVE']=FIVE['ISLET']
open(f"{P}/numbers_extra.tex","w").write("% generated by scripts/37_review_analyses.py\n"+"\n".join(f"\\newcommand{{\\{k}}}{{{v}}}" for k,v in M.items())+"\n")
# ---- results table (with bootstrap support)
L=["\\begin{table}[htbp]\\centering\\footnotesize\\setlength{\\tabcolsep}{4pt}",
"\\caption{Trait--tissue recovery. \\textbf{Rank} is the position of the expected tissue among 11 tissues, ordered by mean enrichment over background. \\textbf{Boot.} is the proportion of 2{,}000 resamples of the trait's variants in which the expected tissue ranks first. $z_{\\text{exp}}$ is the expected tissue's enrichment; $P$ is a one-sided Mann--Whitney test of that tissue against background. Type 2 diabetes appears twice, once per GWAS.}\\label{tab:results}",
"\\begin{tabular}{llrrrlrr}\\toprule","Trait & Expected & $n$ & Rank & Boot. & Top-ranked & $z_{\\text{exp}}$ & $P$ \\\\ \\midrule"]
for t,r in R.sort_values(['rank','visibility'],ascending=[True,False]).iterrows():
    p=f"{r.p_exp:.1e}" if r.p_exp<0.01 else f"{r.p_exp:.2f}"; bold=(lambda s:f"\\textbf{{{s}}}") if t.startswith("Type 2") else (lambda s:s)
    L.append(f"{bold(esc(lab(t)))} & {NICE[r.expected]} & {int(r.n)} & {bold(int(r['rank']))} & {r.boot:.2f} & {NICE[r.top]} & {r.exp_z:.2f} & {p} \\\\")
L+=["\\bottomrule\\end{tabular}\\end{table}"]; open(f"{P}/results_table.tex","w").write("\n".join(L))
# ---- traits table (with GWAS size, cases, fine-mapping regime)
L=["\\begingroup\\footnotesize\\setlength{\\tabcolsep}{4pt}",
"\\begin{longtable}{>{\\raggedright\\arraybackslash}p{3.4cm}>{\\raggedright\\arraybackslash}p{1.4cm}>{\\raggedright\\arraybackslash}p{2.5cm}>{\\raggedright\\arraybackslash}p{1.9cm}r>{\\raggedright\\arraybackslash}p{2.0cm}r}",
"\\caption{Trait panel, data sources and pre-specified expected tissue. $N$ is the GWAS sample size with the number of cases in parentheses for binary traits; fine-mapping regime is as reported by the Open Targets Platform (PICS around reported top hits when no summary statistics were ingested; PICS on summary statistics; SuSiE-inf, with the share of credible sets so mapped where mixed); $n$ is the number of variants scored.}\\label{tab:traits}\\\\",
"\\toprule Trait & Expected & Accession & Source & $N$ (cases) & Fine-mapping & $n$ \\\\ \\midrule \\endfirsthead",
"\\toprule Trait & Expected & Accession & Source & $N$ (cases) & Fine-mapping & $n$ \\\\ \\midrule \\endhead"]
for t,r in R.sort_values('trait').iterrows():
    sid=SID[t]; sm=FM[sid]
    N=f"{sm['nSamples']:,}"+(f" ({sm['nCases']:,})" if sm['nCases'] else "")
    L.append(f"{esc(lab(t))} & {NICE[r.expected]} & \\texttt{{{esc(sid).replace(chr(92)+'_',chr(92)+'_'+chr(92)+'allowbreak ')}}} & {esc(r.source)} & {N} & {reg[t]} & {int(r.n)} \\\\")
L+=["\\bottomrule\\end{longtable}\\endgroup"]; open(f"{P}/traits_table.tex","w").write("\n".join(L))
# ---- track census (with five-mark column)
tc={m:pd.read_csv(f"{D}/tracks_{m}.csv")['grp'].value_counts() for m in ['rna_seq','atac','dnase','chip_histone']}
TC=pd.DataFrame(tc).fillna(0).astype(int).reindex(GR); TC['acc']=TC['atac']+TC['dnase']
L=["\\begin{table}[htbp]\\centering\\small",
"\\caption{AlphaGenome human track census by tissue group. \\textbf{Histone (5 marks)} counts only the five marks scored in this study (H3K27ac, H3K36me3, H3K4me1, H3K4me3, H3K9ac); this is the count used in the track-coverage correlations. The pancreatic islet has no ATAC tracks and a single DNase track.}\\label{tab:tracks}",
"\\begin{tabular}{lrrrrrr}\\toprule","Tissue & RNA-seq & ATAC & DNase & Accessibility & Histone (all) & Histone (5 marks) \\\\ \\midrule"]
for g in GR:
    b=(lambda s:f"\\textbf{{{s}}}") if g=='ISLET' else (lambda s:s)
    L.append(f"{b(NICE[g])} & {b(TC.loc[g,'rna_seq'])} & {b(TC.loc[g,'atac'])} & {b(TC.loc[g,'dnase'])} & {b(TC.loc[g,'acc'])} & {b(TC.loc[g,'chip_histone'])} & {b(FIVE[g])} \\\\")
L+=["\\bottomrule\\end{tabular}\\end{table}"]; open(f"{P}/tracks_table.tex","w").write("\n".join(L))
# ---- relabel LDSC table
s=open(f"{P}/ldsc_table.tex").read()
for k,v in LABEL.items(): s=s.replace(esc(k),esc(v))
open(f"{P}/ldsc_table.tex","w").write(s)
# ---- heatmap with display labels
Mx=pd.read_csv(f"{D}/bench_matrix.csv",index_col=0); ordr=R.sort_values('rank').index.tolist(); Mo=Mx.loc[[t for t in ordr if t in Mx.index]]
fig,ax=plt.subplots(figsize=(7.6,0.34*len(Mo)+1.5)); v=np.nanpercentile(np.abs(Mo.values),97) or 1
im=ax.imshow(Mo.values,cmap='RdBu_r',norm=TwoSlopeNorm(vcenter=0,vmin=-v,vmax=v),aspect='auto')
ax.set_xticks(range(len(GR))); ax.set_xticklabels([NICE[g] for g in GR],rotation=45,ha='right',fontsize=8)
ax.set_yticks(range(len(Mo))); ax.set_yticklabels([lab(t) for t in Mo.index],fontsize=8)
for i,t in enumerate(Mo.index): ax.add_patch(plt.Rectangle((GR.index(R.loc[t,'expected'])-.5,i-.5),1,1,fill=False,ec='black',lw=1.6))
cb=fig.colorbar(im,ax=ax,shrink=.6); cb.set_label('mean enrichment $z$ vs background',fontsize=8); cb.ax.tick_params(labelsize=7)
ax.set_title('Predicted tissue enrichment by trait\n(black box = pre-specified expected tissue)',fontsize=9)
fig.tight_layout(); fig.savefig(f"{P}/fig_heatmap.pdf"); fig.savefig(f"{P}/fig_heatmap.png",dpi=200)
V.drop(columns=['visz']).to_csv(f"{D}/pervariant_summary.csv",index=False); R.to_csv(f"{D}/bench_final_review.csv")
print("\n".join(f"{k} = {v}" for k,v in M.items()))
