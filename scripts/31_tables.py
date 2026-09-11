import os,sys,json,numpy as np,pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
SP=sys.argv[1]; P=f"{SP}/paper"
R=pd.read_csv(f"{SP}/bench_final.csv").rename(columns={"Unnamed: 0":"trait"}); M=pd.read_csv(f"{SP}/bench_matrix.csv",index_col=0)
SM=json.load(open(f"{SP}/study_meta.json"))
GR=list(M.columns)
NICE={'ISLET':'Islet','LIVER':'Liver','ADIPOSE':'Fat','MUSCLE':'Muscle','BRAIN':'Brain',
 'IMMUNE':'Immune','HEART':'Heart','ARTERY':'Artery','LUNG':'Lung','INTESTINE':'Intestine','SKIN':'Skin'}
esc=lambda s:str(s).replace('&','\\&').replace('_','\\_').replace('%','\\%')

# ---- results table
lines=["\\begin{table}[htbp]\\centering\\small",
"\\caption{Trait--tissue recovery. \\textbf{Rank} is the position of the expected tissue among "
+str(len(GR))+" tissues, ordered by mean enrichment over background. $z_{\\text{exp}}$ is the "
"expected tissue's enrichment; $P$ is a one-sided Mann--Whitney test of that tissue against "
"background.}\\label{tab:results}",
"\\begin{tabular}{llrrlrr}\\toprule",
"Trait & Expected & $n$ & Rank & Top-ranked & $z_{\\text{exp}}$ & $P$ \\\\ \\midrule"]
for r in R.sort_values(['rank','trait']).itertuples():
    p=f"{r.p_exp:.1e}" if r.p_exp<0.01 else f"{r.p_exp:.2f}"
    bold=(lambda s:f"\\textbf{{{s}}}") if str(r.trait).startswith("Type 2") else (lambda s:s)
    lines.append(f"{bold(esc(r.trait))} & {NICE[r.expected]} & {r.n} & {bold(r.rank)} & "
                 f"{NICE[r.top]} & {r.exp_z:.2f} & {p} \\\\")
lines+=["\\bottomrule\\end{tabular}\\end{table}"]
open(f"{P}/results_table.tex","w").write("\n".join(lines))

# ---- traits table
ACC=json.load(open(f"{SP}/accessions.json"))
lines=["\\begin{longtable}{p{4.3cm}p{1.5cm}p{3.1cm}rp{2.6cm}}",
"\\caption{Trait panel, data sources and pre-specified expected tissue.}\\label{tab:traits}\\\\",
"\\toprule Trait & Expected & GWAS accession & $N$ & First author \\\\ \\midrule \\endfirsthead",
"\\toprule Trait & Expected & GWAS accession & $N$ & First author \\\\ \\midrule \\endhead"]
for r in R.sort_values('trait').itertuples():
    lines.append(f"{esc(r.trait)} & {NICE[r.expected]} & \\texttt{{\\footnotesize {ACC.get(r.trait,'--')}}} & {r.n} & {esc(r.source)} \\\\")
lines.append("\\bottomrule\\end{longtable}")
open(f"{P}/traits_table.tex","w").write("\n".join(lines))

# ---- track census table
tc={}
for mod in ['rna_seq','atac','dnase','chip_histone']:
    d=pd.read_csv(f"{SP}/tracks_{mod}.csv"); tc[mod]=d['grp'].value_counts()
TC=pd.DataFrame(tc).fillna(0).astype(int).reindex(GR)
TC['accessibility']=TC['atac']+TC['dnase']
lines=["\\begin{table}[htbp]\\centering\\small",
"\\caption{AlphaGenome human track census by tissue group. The pancreatic islet has no ATAC "
"tracks and a single DNase track.}\\label{tab:tracks}",
"\\begin{tabular}{lrrrrr}\\toprule",
"Tissue & RNA-seq & ATAC & DNase & Accessibility & Histone \\\\ \\midrule"]
for g in GR:
    b=(lambda s:f"\\textbf{{{s}}}") if g=='ISLET' else (lambda s:s)
    lines.append(f"{b(NICE[g])} & {b(TC.loc[g,'rna_seq'])} & {b(TC.loc[g,'atac'])} & "
                 f"{b(TC.loc[g,'dnase'])} & {b(TC.loc[g,'accessibility'])} & {b(TC.loc[g,'chip_histone'])} \\\\")
lines+=["\\bottomrule\\end{tabular}\\end{table}"]
open(f"{P}/tracks_table.tex","w").write("\n".join(lines))

# ---- heatmap
ordr=R.sort_values('rank')['trait'].tolist()
Mo=M.loc[[t for t in ordr if t in M.index]]
fig,ax=plt.subplots(figsize=(7.6,0.34*len(Mo)+1.5))
v=np.nanpercentile(np.abs(Mo.values),97) or 1
im=ax.imshow(Mo.values,cmap='RdBu_r',norm=TwoSlopeNorm(vcenter=0,vmin=-v,vmax=v),aspect='auto')
ax.set_xticks(range(len(GR))); ax.set_xticklabels([NICE[g] for g in GR],rotation=45,ha='right',fontsize=8)
ax.set_yticks(range(len(Mo))); ax.set_yticklabels(Mo.index,fontsize=8)
for i,t in enumerate(Mo.index):
    exp=R.set_index('trait').loc[t,'expected']
    ax.add_patch(plt.Rectangle((GR.index(exp)-.5,i-.5),1,1,fill=False,ec='black',lw=1.6))
cb=fig.colorbar(im,ax=ax,shrink=.6); cb.set_label('mean enrichment $z$ vs background',fontsize=8)
cb.ax.tick_params(labelsize=7)
ax.set_title('Predicted tissue enrichment by trait\n(black box = pre-specified expected tissue)',fontsize=9)
fig.tight_layout(); fig.savefig(f"{P}/fig_heatmap.pdf"); fig.savefig(f"{P}/fig_heatmap.png",dpi=200)
print("tables + figure written")
print(TC.to_string())
