import os,time
from alphagenome.models import dna_client, variant_scorers
from alphagenome.data import genome
m=dna_client.create(os.environ['ALPHAGENOME_KEY'])
W=dna_client.SUPPORTED_SEQUENCE_LENGTHS['SEQUENCE_LENGTH_1MB']
v=genome.Variant(chromosome='chr10',position=112998590,reference_bases='C',alternate_bases='T')
iv=v.reference_interval.resize(W)
sc=[variant_scorers.RECOMMENDED_VARIANT_SCORERS[k] for k in ('RNA_SEQ','CHIP_HISTONE')]
t=time.time()
out=m.score_variant(interval=iv,variant=v,variant_scorers=sc)
print(f"elapsed {time.time()-t:.1f}s ; n outputs={len(out)}")
for a in out:
    print("---", type(a).__name__, a.shape)
    print("obs cols:", list(a.obs.columns)[:12])
    print("var cols:", list(a.var.columns)[:12])
    print(a.obs.head(2).to_string())
    print(a.var.head(2).to_string())
