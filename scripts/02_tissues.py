import os, pandas as pd
from alphagenome.models import dna_client
model = dna_client.create(os.environ['ALPHAGENOME_KEY'])
md = model.output_metadata(organism=dna_client.Organism.HOMO_SAPIENS)

df = md.rna_seq
print("rna_seq rows:", len(df))
print("\nunique gtex_tissue (n=%d):" % df['gtex_tissue'].nunique(dropna=True))
print(sorted(df['gtex_tissue'].dropna().unique().tolist()))
print("\nsample biosample_name values (first 40 unique):")
u = sorted(df['biosample_name'].dropna().astype(str).unique().tolist())
print(len(u), "unique biosample_name")
print(u[:40])
print("\ngrep pancrea/islet across biosample_name:")
print([x for x in u if 'pancrea' in x.lower() or 'islet' in x.lower()])
