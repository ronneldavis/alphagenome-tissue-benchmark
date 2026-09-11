import gzip,sys,numpy as np,pandas as pd
L=sys.argv[1]
GRP={
 'ISLET':['Pancreas'],                       # GTEx has no islet; bulk pancreas is the only proxy
 'LIVER':['Liver'],
 'ADIPOSE':['Adipose - Subcutaneous','Adipose - Visceral (Omentum)'],
 'MUSCLE':['Muscle - Skeletal'],
 'BRAIN':None,   # all 'Brain - *'
 'IMMUNE':['Whole Blood','Spleen','Cells - EBV-transformed lymphocytes'],
 'HEART':['Heart - Atrial Appendage','Heart - Left Ventricle'],
 'ARTERY':['Artery - Aorta','Artery - Coronary','Artery - Tibial'],
 'LUNG':['Lung'],
 'INTESTINE':['Colon - Sigmoid','Colon - Transverse','Small Intestine - Terminal Ileum'],
 'SKIN':['Skin - Sun Exposed (Lower leg)','Skin - Not Sun Exposed (Suprapubic)','Cells - Cultured fibroblasts'],
}
df=pd.read_csv(f"{L}/ref/gtex_median_tpm.gct.gz",sep='\t',skiprows=2)
df['gene']=df['Name'].str.split('.').str[0]
tis=[c for c in df.columns if c not in ('Name','Description','gene')]
GRP['BRAIN']=[c for c in tis if c.startswith('Brain - ')]
missing={g:[t for t in ts if t not in tis] for g,ts in GRP.items()}
bad={g:m for g,m in missing.items() if m}
if bad: print("  WARNING unmatched tissue names:",bad)
X=np.log2(df[tis].values+1.0)
coord=pd.read_csv(f"{L}/ref/ENSG_coord_hg19.txt",sep='\t')
ok=df['gene'].isin(set(coord.GENE)).values
X=X[ok]; genes=df['gene'].values[ok]
print(f"  genes with hg19 coords: {len(genes)}  | tissues: {len(tis)}")
# group-level expression = mean of member tissues
G=list(GRP); M=np.column_stack([X[:,[tis.index(t) for t in GRP[g]]].mean(1) for g in G])
# specificity: (group - mean of other groups) / sd across groups  (t-like, Finucane-style)
sd=M.std(1,ddof=1); sd[sd==0]=np.nan
spec=np.column_stack([(M[:,i]-np.delete(M,i,axis=1).mean(1))/sd for i in range(len(G))])
top=int(round(0.10*len(genes)))
print(f"  top 10% per group = {top} genes")
import os; os.makedirs(f"{L}/genesets",exist_ok=True)
for i,g in enumerate(G):
    s=spec[:,i].copy(); s[np.isnan(s)]=-np.inf
    sel=genes[np.argsort(-s)[:top]]
    with open(f"{L}/genesets/{g}.GeneSet","w") as f: f.write("\n".join(sel)+"\n")
    print(f"    {g:10s} n={len(sel)}  tissues={len(GRP[g])}")
# control: all genes
with open(f"{L}/genesets/control.GeneSet","w") as f: f.write("\n".join(genes)+"\n")
print(f"    {'control':10s} n={len(genes)} (all genes, for the control annotation)")
