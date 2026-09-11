import os, pandas as pd
from alphagenome.models import dna_client
pd.set_option('display.width', 200)

model = dna_client.create(os.environ['ALPHAGENOME_KEY'])
md = model.output_metadata(organism=dna_client.Organism.HOMO_SAPIENS)

print("=== output types available ===")
rows=[]
for name in dir(md):
    if name.startswith('_'): continue
    try: v = getattr(md, name)
    except Exception: continue
    if isinstance(v, pd.DataFrame):
        rows.append((name, len(v)))
for n,c in rows: print(f"{n:24s} {c:6d} tracks")

print("\n=== columns on rna_seq ===")
print(list(md.rna_seq.columns))

print("\n=== anything matching pancreas/islet/beta ===")
for n,_ in rows:
    df = getattr(md, n)
    cols = [c for c in df.columns if df[c].dtype==object]
    if not cols: continue
    hay = df[cols].astype(str).agg(' | '.join, axis=1).str.lower()
    hit = df[hay.str.contains('pancrea|islet|langerhans', na=False, regex=True)]
    if len(hit):
        keep=[c for c in ('name','biosample_name','biosample_type','ontology_curie','gtex_tissue','assay') if c in hit.columns]
        print(f"\n--- {n}: {len(hit)} hits ---")
        print(hit[keep].drop_duplicates().to_string(max_rows=40))
