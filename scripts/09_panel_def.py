import sys,pandas as pd
SP=sys.argv[1]
rna=pd.read_csv(f"{SP}/rna_tracks.csv"); hist=pd.read_csv(f"{SP}/hist_tracks.csv")
GROUPS={
 'ISLET'  : r'^type B pancreatic cell$|^endocrine pancreas$|^progenitor cell of endocrine pancreas$',
 'PANCREAS': r'^pancreas$|^body of pancreas$',
 'LIVER'  : r'^liver$|^right lobe of liver$|^hepatocyte$',
 'ADIPOSE': r'adipose|adipocyte|omental fat pad|subcutaneous adipose',
 'MUSCLE' : r'^skeletal muscle|muscle of leg|muscle of trunk|myotube|^myocyte',
}
def tag(df):
    bs=df['biosample_name'].astype(str)
    gt=df['gtex_tissue'].astype(str) if 'gtex_tissue' in df.columns else pd.Series(['']*len(df))
    g=pd.Series([None]*len(df))
    for k,p in GROUPS.items():
        mask=bs.str.contains(p,case=False,regex=True,na=False)
        if k=='LIVER': mask|= gt.eq('Liver')
        if k=='PANCREAS': mask|= gt.eq('Pancreas')
        if k=='ADIPOSE': mask|= gt.str.startswith('Adipose')
        if k=='MUSCLE': mask|= gt.eq('Muscle_Skeletal')
        g[mask.values]=k
    return g
rna['grp']=tag(rna); hist['grp']=tag(hist)
print("=== RNA tracks per group ==="); print(rna['grp'].value_counts().to_string())
print("\n=== HISTONE tracks per group ==="); print(hist['grp'].value_counts().to_string())
print("\n=== histone marks available per group ===")
t=pd.crosstab(hist['grp'],hist['histone_mark'])
common=[c for c in t.columns if (t[c]>0).all()]
print(t[common].to_string() if common else "NO COMMON MARK")
print("\ncommon marks across all 5 groups:",common)
rna.to_csv(f"{SP}/rna_tracks_grp.csv",index=False); hist.to_csv(f"{SP}/hist_tracks_grp.csv",index=False)
print("\nISLET rna tracks:"); print(rna[rna.grp=='ISLET'][['name','biosample_name']].to_string(index=False))
