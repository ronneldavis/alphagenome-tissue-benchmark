import os, re, pandas as pd
from alphagenome.models import dna_client
model = dna_client.create(os.environ['ALPHAGENOME_KEY'])
md = model.output_metadata(organism=dna_client.Organism.HOMO_SAPIENS)

MODS = ['rna_seq','cage','atac','dnase','chip_histone','chip_tf']
PAT = {
 'BETA_CELL'  : r'^type B pancreatic cell$|^endocrine pancreas$',
 'PANCREAS'   : r'^pancreas$|^body of pancreas$|Pancreas',
 'LIVER'      : r'^liver$|^right lobe of liver$|^hepatocyte$|Liver|^HepG2$',
 'ADIPOSE'    : r'adipose|adipocyte|omental fat|subcutaneous',
 'MUSCLE'     : r'skeletal muscle|muscle of leg|muscle of trunk|myotube|myocyte|Muscle_Skeletal',
}
out=[]
for m in MODS:
    df = getattr(md, m)
    bs = df['biosample_name'].astype(str)
    gt = df['gtex_tissue'].astype(str) if 'gtex_tissue' in df.columns else pd.Series(['']*len(df))
    for label,pat in PAT.items():
        hit = df[bs.str.contains(pat,case=False,regex=True,na=False) | gt.str.contains(pat,case=False,regex=True,na=False)]
        out.append({'modality':m,'tissue':label,'n_tracks':len(hit)})
piv = pd.DataFrame(out).pivot(index='tissue',columns='modality',values='n_tracks')
print(piv.to_string())

print("\n=== exact beta-cell tracks per modality ===")
for m in MODS:
    df = getattr(md,m)
    bs = df['biosample_name'].astype(str)
    hit = df[bs.str.contains(r'^type B pancreatic cell$|^endocrine pancreas$|^progenitor cell of endocrine pancreas$',regex=True,na=False)]
    if len(hit):
        print(f"\n{m}: {len(hit)}")
        print(hit[['name','biosample_name','Assay title' if 'Assay title' in hit.columns else 'biosample_type']].drop_duplicates().to_string(max_rows=25,index=False))
