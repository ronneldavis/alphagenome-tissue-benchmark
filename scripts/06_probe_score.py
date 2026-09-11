import os,time
from alphagenome.models import dna_client, variant_scorers
from alphagenome.data import genome
m=dna_client.create(os.environ['ALPHAGENOME_KEY'])
print("=== recommended scorers ===")
for k,v in variant_scorers.RECOMMENDED_VARIANT_SCORERS.items():
    print(f"  {k:28s} {type(v).__name__}  out={getattr(v,'requested_output',None)}")
print("\n=== supported widths ===", dna_client.SUPPORTED_SEQUENCE_LENGTHS)
