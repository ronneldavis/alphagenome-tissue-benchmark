import os,sys,json,gzip,subprocess,time
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
L=sys.argv[1]
B="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
JOBS={
 "GCST90018942":dict(u=f"{B}/GCST90018001-GCST90019000/GCST90018942/harmonised/34594039-GCST90018942-EFO_0004533-Build37.f.tsv.gz",snp="RSID_ADDED",sig="beta,0",n=("N","463178"),addrs=1),
 "GCST90480318":dict(u=f"{B}/GCST90480001-GCST90481000/GCST90480318/harmonised/GCST90480318.h.tsv.gz",snp="rsid",sig="odds_ratio,1",n=("col","n")),
 "GCST90479877":dict(u=f"{B}/GCST90479001-GCST90480000/GCST90479877/GCST90479877.tsv.gz",snp="rsid",sig="odds_ratio,1",n=("col","n")),
 "GCST90480666":dict(u=f"{B}/GCST90480001-GCST90481000/GCST90480666/GCST90480666.tsv.gz",snp="rsid",sig="beta,0",n=("col","n")),
 "GCST90132314":dict(u=f"{B}/GCST90132001-GCST90133000/GCST90132314/harmonised/GCST90132314.h.tsv.gz",snp="variant_id",sig="beta,0",n=("N","1165690")),
 "GCST90000066":dict(u=f"{B}/GCST90000001-GCST90001000/GCST90000066/harmonised/33230300-GCST90000066-EFO_0006335.h.tsv.gz",snp="variant_id",sig="beta,0",n=("col","n")),
 "GCST90244094":dict(u=f"{B}/GCST90244001-GCST90245000/GCST90244094/harmonised/GCST90244094.h.tsv.gz",snp="variant_id",sig="Zscore,0",n=("col","N")),
 "GCST90503108":dict(u=f"{B}/GCST90503001-GCST90504000/GCST90503108/harmonised/GCST90503108.h.tsv.gz",snp="rsid",sig="beta,0",n=("col","n")),
 "GCST006867":dict(u=f"{B}/GCST006001-GCST007000/GCST006867/harmonised/30054458-GCST006867-EFO_0001360-build37.f.tsv.gz",snp="variant_id",sig="beta,0",n=("col","n")),
}
def nsnp(f):
    if not os.path.exists(f): return 0
    with gzip.open(f,'rt') as fh:
        next(fh); return sum(1 for _ in fh)
def dl(a):
    j=JOBS[a]; raw=f"{L}/raw/{a}.tsv.gz"
    for attempt in range(6):
        r=subprocess.run(["curl","-fsL","-C","-","--retry","5","--retry-delay","5",
                          "--connect-timeout","30","--speed-time","120","--speed-limit","20000",
                          "-o",raw,j["u"]],capture_output=True)
        if r.returncode==0: return a,os.path.getsize(raw)
        time.sleep(5)
    return a,-1
print("PHASE A: downloads (3 parallel, resumable)",flush=True)
todo=[a for a in JOBS if nsnp(f"{L}/munged/{a}.sumstats.gz")<200000]
with ThreadPoolExecutor(3) as ex:
    for f in as_completed([ex.submit(dl,a) for a in todo]):
        a,s=f.result(); print(f"  dl {a}: {'FAIL' if s<0 else f'{s/1e6:.0f} MB'}",flush=True)
print("PHASE B: munge",flush=True)
PM=None
for a in todo:
    j=JOBS[a]; raw=f"{L}/raw/{a}.tsv.gz"
    if not os.path.exists(raw): print(f"  {a}: no file"); continue
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
         "--a1","effect_allele","--a2","other_allele","--p","p_value","--signed-sumstats",j["sig"],
         "--merge-alleles",f"{L}/ref/w_hm3.snplist","--out",out]
    cmd += ["--N-col",j["n"][1]] if j["n"][0]=="col" else ["--N",j["n"][1]]
    subprocess.run(cmd,cwd=f"{L}/src3",capture_output=True,text=True)
    n=nsnp(out+".sumstats.gz")
    print(f"  munge {a}: {n} SNPs {'OK' if n>200000 else 'FAIL'}",flush=True)
    for x in (raw,f"{L}/raw/{a}.rs.tsv.gz"):
        if os.path.exists(x): os.remove(x)
print("PHASE C: LDSC-SEG (4 parallel)",flush=True)
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
    return a,"FAILED"
ready=[a for a in json.load(open(f"{L}/sumstats_urls.json")) if nsnp(f"{L}/munged/{a}.sumstats.gz")>200000]
print(f"  traits ready: {len(ready)}",flush=True)
with ThreadPoolExecutor(4) as ex:
    for f in as_completed([ex.submit(cts,a) for a in ready]):
        a,s=f.result(); print(f"  cts {a}: {s}",flush=True)
print("PIPELINE COMPLETE",flush=True)
