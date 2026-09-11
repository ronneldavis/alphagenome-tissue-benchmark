import sys,os,gzip,numpy as np,pandas as pd
L=sys.argv[1]
G=['ISLET','LIVER','ADIPOSE','MUSCLE','BRAIN','IMMUNE','HEART','ARTERY','LUNG','INTESTINE','SKIN']
OUT=f"{L}/cts_ldscores"; os.makedirs(OUT,exist_ok=True)
for c in range(1,23):
    df=pd.read_csv(f"{L}/annot/cts.{c}.l2.ldscore.gz",sep='\t')
    M=open(f"{L}/annot/cts.{c}.l2.M").read().split()
    M5=open(f"{L}/annot/cts.{c}.l2.M_5_50").read().split()
    names=['base','control']+G                      # annot column order
    assert len(M5)==len(names), (len(M5),len(names))
    for i,nm in enumerate(names):
        if nm=='base': continue
        sub=df[['CHR','SNP','BP',nm+'L2']].copy()
        sub.to_csv(f"{OUT}/{nm}.{c}.l2.ldscore.gz",sep='\t',index=False,compression='gzip')
        open(f"{OUT}/{nm}.{c}.l2.M","w").write(M[i]+"\n")
        open(f"{OUT}/{nm}.{c}.l2.M_5_50","w").write(M5[i]+"\n")
    if c in (1,22): print(f"  chr{c}: split {len(names)-1} annotations, {len(df)} SNPs")
with open(f"{L}/tissues.ldcts","w") as f:
    for g in G: f.write(f"{g}\t{OUT}/{g}.,{OUT}/control.\n")
print("  ldcts file:"); print(open(f"{L}/tissues.ldcts").read().strip().split("\n")[0][:120])
