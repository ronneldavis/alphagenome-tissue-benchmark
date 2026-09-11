import os,sys,pandas as pd
from alphagenome.models import dna_client
md=dna_client.create(os.environ['ALPHAGENOME_KEY']).output_metadata(organism=dna_client.Organism.HOMO_SAPIENS)
G={
 'ISLET'  : (r'^type B pancreatic cell$|^endocrine pancreas$|^progenitor cell of endocrine pancreas$', None),
 'LIVER'  : (r'^liver$|^right lobe of liver$|^hepatocyte$|^HepG2$', 'Liver'),
 'ADIPOSE': (r'adipose|adipocyte|omental fat pad', 'Adipose'),
 'MUSCLE' : (r'^skeletal muscle|muscle of leg|muscle of trunk|myotube|^myocyte', 'Muscle_Skeletal'),
 'BRAIN'  : (r'brain|cortex|neuron|astrocyte|hippocamp|cerebell|caudate|putamen|substantia nigra|spinal cord', 'Brain'),
 'IMMUNE' : (r'T cell|B cell|monocyte|macrophage|natural killer|lymphocyte|CD4|CD8|neutrophil|GM12878|dendritic|thymus|spleen', 'Whole_Blood'),
 'HEART'  : (r'heart|cardiac|ventricle|atrium|cardiomyocyte', 'Heart'),
 'ARTERY' : (r'artery|aorta|endothelial cell of|smooth muscle cell|vascular', 'Artery'),
 'LUNG'   : (r'^lung$|bronchial|alveolar|^upper lobe of left lung$|A549', 'Lung'),
 'INTESTINE':(r'intestine|colon|sigmoid|ileum|HT-29|Caco-2', 'Colon'),
 'SKIN'   : (r'^skin|keratinocyte|foreskin|fibroblast of dermis|melanocyte', 'Skin'),
}
def tag(df):
    bs=df['biosample_name'].astype(str)
    gt=df['gtex_tissue'].astype(str) if 'gtex_tissue' in df.columns else pd.Series(['']*len(df))
    out=pd.Series([None]*len(df))
    for k,(pat,gpre) in G.items():
        m=bs.str.contains(pat,case=False,regex=True,na=False)
        if gpre: m|=gt.str.startswith(gpre)
        out[m.values]=k
    return out
rows={}
for mod in ['rna_seq','atac','dnase','chip_histone','cage']:
    df=getattr(md,mod).copy(); df['grp']=tag(df); rows[mod]=df['grp'].value_counts()
    df.to_csv(f"{sys.argv[1]}/tracks_{mod}.csv",index=False)
T=pd.DataFrame(rows).fillna(0).astype(int)
print(T.to_string())
h=pd.read_csv(f"{sys.argv[1]}/tracks_chip_histone.csv")
ct=pd.crosstab(h['grp'],h['histone_mark'])
common=[c for c in ct.columns if (ct[c]>0).all()]
print("\nhistone marks in ALL groups:",common)
print("\ngroups missing RNA or histone:",[g for g in G if T.loc[g,'rna_seq']==0 or T.loc[g,'chip_histone']==0])
