import os,sys,gzip,subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
L=sys.argv[1]
B="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
# Derived from 90b_pipeline_v2.py. Downloads for these 3 accessions are already
# complete and verified in raw/ (size == remote Content-Length, gzip -t OK).
# Same column settings as 90b_pipeline_v2.py JOBS. Runs PHASE B (munge) then
# PHASE C (--h2-cts) only -- no PHASE A download step.
JOBS={
 "GCST90446794":dict(  # Ulcerative colitis (Liu 2023, Nat Genet) -> IMMUNE
   u=f"{B}/GCST90446001-GCST90447000/GCST90446794/harmonised/GCST90446794.h.tsv.gz",
   snp="rsid",a1="effect_allele",a2="other_allele",p="p_value",
   sig="beta,0",n=("N","375508")),
 "GCST90013445":dict(  # Type 1 diabetes (Robertson 2021, Nat Genet) -> IMMUNE
   u=f"{B}/GCST90013001-GCST90014000/GCST90013445/harmonised/GCST90013445.h.tsv.gz",
   snp="rsid",a1="effect_allele",a2="other_allele",p="p_value",
   sig="beta,0",n=("N","59527")),
 "GCST006250":dict(  # Intelligence (Savage 2018, Nat Genet) -> BRAIN
   u=f"{B}/GCST006001-GCST007000/GCST006250/harmonised/29942086-GCST006250-EFO_0004337.h.tsv.gz",
   snp="hm_rsid",a1="hm_effect_allele",a2="hm_other_allele",p="p_value",
   sig="hm_beta,0",n=("col","n_analyzed"),
   # FIX (2026-09-10): this harmonised file also carries the original unprefixed
   # effect_allele/other_allele columns. munge_sumstats.py's default cname
   # dictionary auto-maps those to A1/A2 IN ADDITION to our explicit --a1/--a2
   # (hm_effect_allele/hm_other_allele), producing duplicate A1/A2 columns and
   # "AttributeError: 'DataFrame' object has no attribute 'str'". --ignore
   # drops the unprefixed duplicates before column matching.
   ignore="effect_allele,other_allele"),
}
# Optional: restrict this run to a single accession (used for column-name fix
# reruns) so it never touches another accession's in-flight munge/cts.
ACC_FILTER=sys.argv[2] if len(sys.argv)>2 else None
def nsnp(f):
    if not os.path.exists(f): return 0
    with gzip.open(f,'rt') as fh:
        next(fh); return sum(1 for _ in fh)
print("PHASE B: munge (downloads already complete, skipping PHASE A)",flush=True)
PM=None
todo=[ACC_FILTER] if ACC_FILTER else list(JOBS.keys())
for a in todo:
    j=JOBS[a]; raw=f"{L}/raw/{a}.tsv.gz"
    if nsnp(f"{L}/munged/{a}.sumstats.gz")>200000:
        print(f"  {a}: already munged, skip",flush=True); continue
    if not os.path.exists(raw): print(f"  {a}: no raw file",flush=True); continue
    src=raw
    if j.get("addrs"):
        if PM is None:
            fr=[pd.read_csv(f"{L}/ref/1000G_EUR_Phase3_plink/1000G.EUR.QC.{c}.bim",sep='\t',header=None,
                 names=['CHR','SNP','CM','BP','A1','A2'],usecols=['CHR','SNP','BP']) for c in range(1,23)]
            m=pd.concat(fr); m['k']=m.CHR.astype(str)+":"+m.BP.astype(str); PM=dict(zip(m.k,m.SNP))
        src=f"{L}/raw/{a}.rs.tsv.gz"; first=True
        with gzip.open(src,"wt") as o:
            for ch in pd.read_csv(raw,sep="\t",chunksize=1_000_000,low_memory=False,na_values=[".","NA","#NA",""]):
                k=ch["chromosome"].astype(str)+":"+ch["base_pair_location"].astype("Int64").astype(str)
                ch["RSID_ADDED"]=k.map(PM); ch=ch.dropna(subset=["RSID_ADDED"])
                ch.to_csv(o,sep="\t",index=False,header=first); first=False
    out=f"{L}/munged/{a}"
    cmd=[f"{L}/venv/bin/python","-W","ignore","munge_sumstats.py","--sumstats",src,"--snp",j["snp"],
         "--a1",j["a1"],"--a2",j["a2"],"--p",j["p"],"--signed-sumstats",j["sig"],
         "--merge-alleles",f"{L}/ref/w_hm3.snplist","--out",out]
    cmd += ["--N-col",j["n"][1]] if j["n"][0]=="col" else ["--N",j["n"][1]]
    if j.get("ignore"): cmd += ["--ignore",j["ignore"]]
    r=subprocess.run(cmd,cwd=f"{L}/src3",capture_output=True,text=True)
    n=nsnp(out+".sumstats.gz")
    print(f"  munge {a}: {n} SNPs {'OK' if n>200000 else 'FAIL'}",flush=True)
    if n<=200000:
        print(f"  --- munge {a} stderr tail ---",flush=True)
        print("\n".join(r.stderr.splitlines()[-20:]),flush=True)
    for x in (raw,f"{L}/raw/{a}.rs.tsv.gz"):
        if os.path.exists(x): os.remove(x)
print("PHASE C: LDSC-SEG (3 parallel)",flush=True)
def cts(a):
    o=f"{L}/results/{a}"
    if os.path.exists(o+".cell_type_results.txt"): return a,"cached"
    r=subprocess.run([f"{L}/venv/bin/python","-W","ignore","ldsc.py","--h2-cts",
        f"{L}/munged/{a}.sumstats.gz","--ref-ld-chr",f"{L}/ref/baselineLD.",
        "--ref-ld-chr-cts",f"{L}/tissues.ldcts",
        "--w-ld-chr",f"{L}/ref/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.",
        "--out",o],cwd=f"{L}/src3",capture_output=True,text=True)
    if os.path.exists(o+".cell_type_results.txt"):
        top=open(o+".cell_type_results.txt").readlines()[1].split()
        return a,f"top={top[0]} P={float(top[3]):.3g}"
    tail=r.stderr.splitlines()[-8:] if r.stderr else []
    return a,f"FAILED rc={r.returncode} stderr_tail={tail}"
ready=[a for a in JOBS if nsnp(f"{L}/munged/{a}.sumstats.gz")>200000]
if ACC_FILTER: ready=[a for a in ready if a==ACC_FILTER]
print(f"  traits ready: {len(ready)} -> {ready}",flush=True)
with ThreadPoolExecutor(3) as ex:
    for f in as_completed([ex.submit(cts,a) for a in ready]):
        a,s=f.result(); print(f"  cts {a}: {s}",flush=True)
print("PIPELINE COMPLETE",flush=True)
