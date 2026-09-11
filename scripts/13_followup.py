import sys,numpy as np,pandas as pd
from scipy.stats import mannwhitneyu
SP=sys.argv[1]
GRPS=['ISLET','PANCREAS','LIVER','ADIPOSE','MUSCLE']
T=pd.read_csv(f"{SP}/t2d_tissue_z.csv"); L=pd.read_csv(f"{SP}/ldl_tissue_z.csv")

print("="*78); print("Q1. Which LDL loci drive the liver signal? (top 12)"); print("="*78)
print(L.sort_values('LIVER',ascending=False).head(12)[['rsid','l2g_gene','gene_used','pip','LIVER','ISLET','MUSCLE']].round(2).to_string(index=False))
print("\nmean LIVER z =",round(L.LIVER.mean(),3)," median =",round(L.LIVER.median(),3),
      "| share of loci with z>1:",f"{(L.LIVER>1).mean()*100:.1f}%")

print("\n"+"="*78); print("Q2. T2D: top loci by ISLET score, and where is TCF7L2/rs7903146?"); print("="*78)
print(T.sort_values('ISLET',ascending=False).head(12)[['rsid','l2g_gene','gene_used','pip','ISLET','LIVER','ADIPOSE','MUSCLE']].round(2).to_string(index=False))
T['islet_rank']=T['ISLET'].rank(ascending=False).astype(int)
tc=T[T.l2g_gene=='TCF7L2']
print("\nTCF7L2 signals (the textbook islet locus):")
print(tc[['rsid','pip','islet_rank','ISLET','LIVER','ADIPOSE','MUSCLE']].round(2).to_string(index=False))

print("\n"+"="*78); print("Q3. Do T2D variants show ANY effect above background, in any track?"); print("="*78)
def glob(tag):
    z=np.load(f"{SP}/{tag}_scores.npz")
    R,H=np.abs(z['rna']),np.abs(z['hist'])
    return R.mean(1), H.mean(1), np.concatenate([R,H],axis=1).max(1)
gb=glob('bg'); gt=glob('t2d'); gl=glob('ldl')
for nm,(a,b,c),lab in ((0,gt,'T2D'),(1,gl,'LDL')):
    pass
for lab,g in (('T2D',gt),('LDL',gl)):
    print(f"\n{lab} vs background:")
    for i,mname in enumerate(['mean|RNA LFC| (all 371 tracks)','mean|histone| (all 1116 tracks)','max|effect| (any track)']):
        p=mannwhitneyu(g[i],gb[i],alternative='greater').pvalue
        print(f"   {mname:36s} trait={g[i].mean():.4f}  bg={gb[i].mean():.4f}  ratio={g[i].mean()/gb[i].mean():.2f}x  p={p:.2e}")
