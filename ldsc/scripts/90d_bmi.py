import os,sys,gzip,subprocess
L=sys.argv[1]
B="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
# Derived from 90b_pipeline_v2.py / 90c_finish.py, for GCST90474606 (BMI) only.
# Its raw/GCST90474606.tsv.gz was finished via a segmented parallel download
# (94_bmi_segments.py) after the single-connection download stalled at
# ~190 KB/s; final size verified == remote Content-Length and gzip -t OK.
# Same column settings as 90b_pipeline_v2.py JOBS. 1 job, munge then h2-cts.
a="GCST90474606"
j=dict(  # Body mass index (UKB WGS Consortium 2025, Nature) -> BRAIN
   u=f"{B}/GCST90474001-GCST90475000/GCST90474606/harmonised/GCST90474606.h.tsv.gz",
   snp="rsid",a1="effect_allele",a2="other_allele",p="p_value",
   sig="beta,0",n=("col","n"))
def nsnp(f):
    if not os.path.exists(f): return 0
    with gzip.open(f,'rt') as fh:
        next(fh); return sum(1 for _ in fh)
raw=f"{L}/raw/{a}.tsv.gz"
out=f"{L}/munged/{a}"
if nsnp(out+".sumstats.gz")>200000:
    print(f"  {a}: already munged, skip",flush=True)
else:
    print("PHASE B: munge",flush=True)
    if not os.path.exists(raw):
        print(f"  {a}: no raw file",flush=True); sys.exit(1)
    cmd=[f"{L}/venv/bin/python","-W","ignore","munge_sumstats.py","--sumstats",raw,"--snp",j["snp"],
         "--a1",j["a1"],"--a2",j["a2"],"--p",j["p"],"--signed-sumstats",j["sig"],
         "--merge-alleles",f"{L}/ref/w_hm3.snplist","--out",out]
    cmd += ["--N-col",j["n"][1]] if j["n"][0]=="col" else ["--N",j["n"][1]]
    r=subprocess.run(cmd,cwd=f"{L}/src3",capture_output=True,text=True)
    n=nsnp(out+".sumstats.gz")
    print(f"  munge {a}: {n} SNPs {'OK' if n>200000 else 'FAIL'}",flush=True)
    if n<=200000:
        print(f"  --- munge {a} stderr tail ---",flush=True)
        print("\n".join(r.stderr.splitlines()[-20:]),flush=True)
        sys.exit(1)
    if os.path.exists(raw): os.remove(raw)
print("PHASE C: LDSC-SEG (1 job)",flush=True)
o=f"{L}/results/{a}"
if os.path.exists(o+".cell_type_results.txt"):
    print(f"  cts {a}: cached",flush=True)
else:
    r=subprocess.run([f"{L}/venv/bin/python","-W","ignore","ldsc.py","--h2-cts",
        f"{L}/munged/{a}.sumstats.gz","--ref-ld-chr",f"{L}/ref/baselineLD.",
        "--ref-ld-chr-cts",f"{L}/tissues.ldcts",
        "--w-ld-chr",f"{L}/ref/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.",
        "--out",o],cwd=f"{L}/src3",capture_output=True,text=True)
    if os.path.exists(o+".cell_type_results.txt"):
        top=open(o+".cell_type_results.txt").readlines()[1].split()
        print(f"  cts {a}: top={top[0]} P={float(top[3]):.3g}",flush=True)
    else:
        tail=r.stderr.splitlines()[-8:] if r.stderr else []
        print(f"  cts {a}: FAILED rc={r.returncode} stderr_tail={tail}",flush=True)
print("PIPELINE COMPLETE",flush=True)
