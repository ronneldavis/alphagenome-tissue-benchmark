import os,sys,time,json
import numpy as np, pandas as pd
from alphagenome.models import dna_client, variant_scorers
from alphagenome.data import genome

SP, TAG = sys.argv[1], sys.argv[2]
df = pd.read_csv(f"{SP}/{TAG}_variants.csv")
df = df[~df.coding].reset_index(drop=True)
print(f"{TAG}: scoring {len(df)} non-coding signals", flush=True)

m = dna_client.create(os.environ['ALPHAGENOME_KEY'])
W = dna_client.SUPPORTED_SEQUENCE_LENGTHS['SEQUENCE_LENGTH_1MB']
SC = [variant_scorers.RECOMMENDED_VARIANT_SCORERS[k] for k in ('RNA_SEQ','CHIP_HISTONE')]

rna_rows, hist_rows, meta, rna_var, hist_var = [], [], [], None, None
t0=time.time()
for i, r in df.iterrows():
    v = genome.Variant(chromosome=f"chr{r.chrom}", position=int(r.pos),
                       reference_bases=str(r.ref), alternate_bases=str(r.alt))
    for attempt in range(4):
        try:
            out = m.score_variant(interval=v.reference_interval.resize(W), variant=v, variant_scorers=SC)
            break
        except Exception as e:
            if attempt==3:
                print(f"  FAIL {r.variant_id}: {type(e).__name__} {str(e)[:100]}", flush=True); out=None; break
            time.sleep(2*(attempt+1))
    if out is None: continue
    rna, hist = out[0], out[1]
    if rna_var is None: rna_var, hist_var = rna.var.copy(), hist.var.copy()
    X = np.asarray(rna.X)                       # genes x tracks
    genes = rna.obs['gene_name'].astype(str).tolist()
    if r.l2g_gene in genes: gi = genes.index(r.l2g_gene)
    else: gi = int(np.nanargmax(np.nanmean(np.abs(X),axis=1))) if X.size else 0
    rna_rows.append(X[gi]); hist_rows.append(np.asarray(hist.X)[0])
    meta.append(dict(variant_id=r.variant_id, rsid=r.rsid, pip=r.pip, cs_size=r.cs_size,
                     l2g_gene=r.l2g_gene, gene_used=genes[gi] if genes else None))
    if (i+1)%25==0: print(f"  {i+1}/{len(df)}  {time.time()-t0:.0f}s", flush=True)

np.savez_compressed(f"{SP}/{TAG}_scores.npz", rna=np.array(rna_rows), hist=np.array(hist_rows))
pd.DataFrame(meta).to_csv(f"{SP}/{TAG}_scored_meta.csv", index=False)
rna_var.to_csv(f"{SP}/rna_tracks.csv", index=False); hist_var.to_csv(f"{SP}/hist_tracks.csv", index=False)
print(f"{TAG}: done {len(rna_rows)} scored in {time.time()-t0:.0f}s", flush=True)
