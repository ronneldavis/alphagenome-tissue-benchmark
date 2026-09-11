import sys,os,gzip,numpy as np,pandas as pd
L=sys.argv[1]; CHR=int(sys.argv[2]); WIN=100_000
G=['ISLET','LIVER','ADIPOSE','MUSCLE','BRAIN','IMMUNE','HEART','ARTERY','LUNG','INTESTINE','SKIN']
coord=pd.read_csv(f"{L}/ref/ENSG_coord_hg19.txt",sep='\t')
coord=coord[coord.CHR.astype(str)==str(CHR)]
bim=pd.read_csv(f"{L}/ref/1000G_EUR_Phase3_plink/1000G.EUR.QC.{CHR}.bim",sep='\t',header=None,
                names=['CHR','SNP','CM','BP','A1','A2'])
bp=bim.BP.values
def annot_for(gs):
    s=set(open(f"{L}/genesets/{gs}.GeneSet").read().split())
    sub=coord[coord.GENE.isin(s)]
    if len(sub)==0: return np.zeros(len(bp),dtype=np.int8)
    iv=sorted(zip((sub.START-WIN).clip(lower=0).values,(sub.END+WIN).values))
    merged=[];cs,ce=iv[0]
    for a,b in iv[1:]:
        if a<=ce: ce=max(ce,b)
        else: merged.append((cs,ce)); cs,ce=a,b
    merged.append((cs,ce))
    starts=np.array([m[0] for m in merged]); ends=np.array([m[1] for m in merged])
    idx=np.searchsorted(starts,bp,side='right')-1
    out=np.zeros(len(bp),dtype=np.int8)
    ok=idx>=0
    out[ok]=(bp[ok]<=ends[idx[ok]]).astype(np.int8)
    return out
cols={'CHR':bim.CHR.values,'BP':bim.BP.values,'SNP':bim.SNP.values,'CM':bim.CM.values,
      'base':np.ones(len(bp),dtype=np.int8),'control':annot_for('control')}
for g in G: cols[g]=annot_for(g)
df=pd.DataFrame(cols)
os.makedirs(f"{L}/annot",exist_ok=True)
df.to_csv(f"{L}/annot/cts.{CHR}.annot.gz",sep='\t',index=False,compression='gzip')
print(f"  chr{CHR}: {len(bp)} SNPs | " + " ".join(f"{g}={int(df[g].sum())}" for g in ['ISLET','LIVER','BRAIN','IMMUNE']) + f" control={int(df['control'].sum())}")
