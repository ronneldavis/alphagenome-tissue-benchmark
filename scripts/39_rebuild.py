"""Rebuild the benchmark + LDSC-comparison outputs from a panel definition file.

This replaces the old hard-coded chain that used to run as

    30_bench.py -> (an un-tracked BH-correction step, scripts/33_t2d_deep.py)
    -> 34_explain.py -> 35_visibility.py -> 36_numbers2.py
    -> ldsc/scripts/95_compare.py -> (a lost "NMATCH-family" script)

against a scratch directory `SP` that held both the score files and the data
CSVs together. Now the score files live in `<P>/scores/` and the data CSVs in
`<P>/data/`, and the trait/study/tissue/source/file mapping that used to be
hard-coded as the `LABEL` dict in `30_bench.py` and the `FILES` dict in
`35_visibility.py` lives in a panel JSON file instead (see `data/panel_v1.json`).

Every formula below is copied unchanged from the original scripts; only the
*source* of the trait/study/tissue/file mapping changed (panel file instead of
hard-coded dict), and the LDSC matched-set mapping now comes from the panel
file instead of the old `ldsc/sumstats_urls.json` (an entry counts as
"matched" iff `<P>/ldsc/results/<study>.cell_type_results.txt` exists).

One exception, called out explicitly because it is a deliberate correction,
not a re-derivation: the 'visibility' column always uses the pooled-track
formula from `35_visibility.py`, never the `bench_matrix.max(axis=1)` version
that `34_explain.py` computed for its own (unused downstream) report -- see
the review notes in REVIEW.md section 1 for why that version was mechanical.

Usage:
    uv run --with numpy --with pandas --with scipy --with statsmodels --with matplotlib \\
        python scripts/39_rebuild.py data/panel_v1.json .

Writes (all under <P>):
    data/bench_summary.csv, data/bench_matrix.csv, data/bench_explain.csv,
    data/bench_final.csv, data/method_comparison.csv, ldsc_table.tex, numbers.tex
"""
import sys, os, json
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, binomtest, false_discovery_control, spearmanr

PANEL_PATH, P = sys.argv[1], sys.argv[2]
SC = f"{P}/scores"
D = f"{P}/data"

MARKS = ['H3K27ac', 'H3K36me3', 'H3K4me1', 'H3K4me3', 'H3K9ac']
GR = ['ISLET', 'LIVER', 'ADIPOSE', 'MUSCLE', 'BRAIN', 'IMMUNE', 'HEART', 'ARTERY', 'LUNG', 'INTESTINE', 'SKIN']
NT = 11
NICE = {'ISLET': 'Islet', 'LIVER': 'Liver', 'ADIPOSE': 'Fat', 'MUSCLE': 'Muscle', 'BRAIN': 'Brain',
        'IMMUNE': 'Immune', 'HEART': 'Heart', 'ARTERY': 'Artery', 'LUNG': 'Lung', 'INTESTINE': 'Intestine', 'SKIN': 'Skin'}
esc = lambda s: str(s).replace('&', '\\&').replace('_', '\\_').replace('%', '\\%')


def sci(p):
    """Same P-value formatter as 36_numbers2.py / 37_review_analyses.py."""
    if p >= 0.01:
        return f"{p:.2f}"
    m, e = f"{p:.1e}".split("e")
    return f"{m}\\times 10^{{{int(e)}}}"


panel = json.load(open(PANEL_PATH))

# =========================================================== 30_bench.py ===
h = pd.read_csv(f"{D}/tracks_chip_histone.csv")
ok = h['histone_mark'].isin(MARKS).values
IDX = {g: ((h['grp'] == g).values & ok) for g in GR}
BG = np.abs(np.load(f"{SC}/bg_scores.npz")['hist'])
NR = {g: (BG[:, IDX[g]].mean(1).mean(), BG[:, IDX[g]].mean(1).std()) for g in GR}

rows = []
mat = {}
for e in panel:
    trait, exp, src, fpfx = e['trait'], e['expected'], e['source'], e['file']
    fpath = f"{SC}/{fpfx}_scores.npz"
    if not os.path.exists(fpath):
        continue
    H = np.abs(np.load(fpath)['hist'])
    if H.ndim != 2 or H.shape[0] < 20:
        continue
    z = {g: (H[:, IDX[g]].mean(1) - NR[g][0]) / NR[g][1] for g in GR}
    means = {g: z[g].mean() for g in GR}
    order = sorted(GR, key=lambda g: -means[g])
    rank = order.index(exp) + 1
    mat[trait] = means
    rows.append(dict(
        trait=trait, source=src, n=H.shape[0], expected=exp, rank=rank,
        top=order[0], top_z=round(means[order[0]], 3), exp_z=round(means[exp], 3),
        p_exp=mannwhitneyu(z[exp], (BG[:, IDX[exp]].mean(1) - NR[exp][0]) / NR[exp][1],
                            alternative='greater').pvalue))
R = pd.DataFrame(rows).sort_values('rank')

# ============================================ BH q-values (33_t2d_deep.py) ==
R['q_exp'] = false_discovery_control(R['p_exp'].values, method='bh')
R.to_csv(f"{D}/bench_summary.csv", index=False)
pd.DataFrame(mat).T[GR].to_csv(f"{D}/bench_matrix.csv")

# ============================ 34_explain.py (n_tracks_exp / acc_tracks_exp) =
# NOTE: 34_explain.py's own 'visibility' (bench_matrix.max(axis=1)) is
# deliberately NOT reproduced here -- see module docstring. Only the track
# census columns are reused, exactly as 35_visibility.py itself reuses them.
NTR = {g: int(IDX[g].sum()) for g in GR}
acc = (pd.read_csv(f"{D}/tracks_atac.csv")['grp'].value_counts()
       .add(pd.read_csv(f"{D}/tracks_dnase.csv")['grp'].value_counts(), fill_value=0)
       .fillna(0))
R2 = R.set_index('trait').copy()
R2['n_tracks_exp'] = [NTR[g] for g in R2['expected']]
R2['acc_tracks_exp'] = [int(acc.get(g, 0) or 0) for g in R2['expected']]
R2['hit1'] = (R2['rank'] == 1).astype(int)
R2.to_csv(f"{D}/bench_explain.csv")

# =================================================== 35_visibility.py =======
ALL = np.zeros(len(h), bool)
for g in GR:
    ALL |= IDX[g]
bgv = BG[:, ALL].mean(1)
mu, sd = bgv.mean(), bgv.std()
FILES = {e['trait']: e['file'] for e in panel}
Rf = pd.read_csv(f"{D}/bench_summary.csv").set_index('trait')
vis = {}
for t, f in FILES.items():
    fpath = f"{SC}/{f}_scores.npz"
    if not os.path.exists(fpath):
        continue
    H = np.abs(np.load(fpath)['hist'])
    vis[t] = ((H[:, ALL].mean(1) - mu) / sd).mean()
Rf['visibility'] = pd.Series(vis)
Rf['recovered'] = (Rf['rank'] == 1).astype(int)
E = pd.read_csv(f"{D}/bench_explain.csv", index_col=0)
Rf['n_tracks_exp'] = E['n_tracks_exp']
Rf['acc_tracks_exp'] = E['acc_tracks_exp']
Rf.to_csv(f"{D}/bench_final.csv")

# =================================================== 36_numbers2.py =========
Rn = pd.read_csv(f"{D}/bench_final.csv").rename(columns={"Unnamed: 0": "trait"})
tc = {m: pd.read_csv(f"{D}/tracks_{m}.csv")['grp'].value_counts() for m in ['rna_seq', 'atac', 'dnase', 'chip_histone']}
TC = pd.DataFrame(tc).fillna(0).astype(int)
nbg = int(np.load(f"{SC}/bg_scores.npz")['hist'].shape[0])
n = len(Rn)
h1 = int((Rn['rank'] == 1).sum())
h3 = int((Rn['rank'] <= 3).sum())
nsig = int((Rn.q_exp < 0.05).sum())
rho_v, p_v = spearmanr(Rn['visibility'], Rn['exp_z'])
rho_t, p_t = spearmanr(Rn['n_tracks_exp'], Rn['exp_z'])
rho_a, p_a = spearmanr(Rn['acc_tracks_exp'], Rn['exp_z'])
a = Rn[Rn.q_exp < 0.05]['visibility']
b = Rn[Rn.q_exp >= 0.05]['visibility']
pmw = mannwhitneyu(a, b, alternative='greater').pvalue
g = Rn.set_index('trait')

M = {
    'NTRAITS': n, 'NBG': nbg, 'NTISSUE': NT, 'NHITONE': h1, 'NHITTHREE': h3, 'NSIG': nsig, 'NCAP': 150,
    'HITONEPCT': f"{h1 / n * 100:.0f}", 'HITTHREEPCT': f"{h3 / n * 100:.0f}",
    'PHITONE': sci(binomtest(h1, n, 1 / NT, alternative='greater').pvalue),
    'PHITTHREE': sci(binomtest(h3, n, 3 / NT, alternative='greater').pvalue),
    'CHANCEONE': f"{n / NT:.1f}", 'CHANCETHREE': f"{3 * n / NT:.1f}",
    'RHOVIS': f"{rho_v:+.2f}", 'PVIS': sci(p_v), 'RHOTRACK': f"{rho_t:+.2f}", 'PTRACK': sci(p_t),
    'RHOACC': f"{rho_a:+.2f}", 'PACC': sci(p_a),
    'VISREC': f"{a.median():.2f}", 'VISNOT': f"{b.median():.2f}", 'PVISMWU': sci(pmw),
    'BRAINHIST': int(TC.loc['BRAIN', 'chip_histone']),
    'BRAINACC': int(TC.loc['BRAIN', 'atac'] + TC.loc['BRAIN', 'dnase']),
    'ISLETRNA': int(TC.loc['ISLET', 'rna_seq']), 'ISLETHIST': int(TC.loc['ISLET', 'chip_histone']),
    'ISLETACC': int(TC.loc['ISLET', 'atac'] + TC.loc['ISLET', 'dnase']),
    'LIVERACC': int(TC.loc['LIVER', 'atac'] + TC.loc['LIVER', 'dnase']),
    'IMMUNEACC': int(TC.loc['IMMUNE', 'atac'] + TC.loc['IMMUNE', 'dnase']),
    'CODEURL': "the repository listed under Data availability",
}
# These macros are about specific named traits from the original 20-pair
# panel (they feed prose in main.tex that names those traits explicitly);
# they are only meaningful while the panel still contains those trait names.
COGT = 'Intelligence' if 'Intelligence' in g.index else 'Cognitive ability / education'  # cognition trait present in this panel
for key, trait in [('TDZ', 'Type 2 diabetes'), ('TDXZ', 'Type 2 diabetes (replication)'),
                    ('BMIZ', 'Body mass index'), ('COGZ', COGT),
                    ('ALTZ', 'Alanine aminotransferase')]:
    if trait in g.index:
        M[key] = f"{g.loc[trait, 'exp_z']:.2f}"
if COGT in g.index:
    M['COGRANK'] = int(g.loc[COGT, 'rank']); M['COGP'] = f"{g.loc[COGT, 'p_exp']:.2f}"
if 'Body mass index' in g.index:
    M['BMIP'] = f"{g.loc['Body mass index', 'p_exp']:.2f}"
if 'Type 2 diabetes' in g.index:
    M['TDP'] = f"{g.loc['Type 2 diabetes', 'p_exp']:.3f}"
if 'Type 2 diabetes (replication)' in g.index:
    M['TDXRANK'] = int(g.loc['Type 2 diabetes (replication)', 'rank'])
    M['TDXP'] = f"{g.loc['Type 2 diabetes (replication)', 'p_exp']:.2f}"
if 'Body mass index' in g.index:
    M['BMIRANK'] = int(g.loc['Body mass index', 'rank'])

# OVERLAP/OVERLAPPCT/CSMEDIAN: hard-coded to the two T2D score files exactly
# as 36_numbers2.py hard-coded them (this is the Mahajan/Xue comparison, not
# a function of the panel).
vm = set(pd.read_csv(f"{SC}/T_GCST009379_variants.csv").variant_id)
vx = set(pd.read_csv(f"{SC}/t2d_variants.csv").variant_id)
M['OVERLAP'] = len(vm & vx)
M['OVERLAPPCT'] = f"{len(vm & vx) / min(len(vm), len(vx)) * 100:.0f}"
M['CSMEDIAN'] = int(pd.read_csv(f"{SC}/t2d_variants.csv")['cs_size'].median())

# ================================== 95_compare.py, adapted to the panel file
# An entry is "matched" (run through both AlphaGenome and LDSC-SEG) iff its
# study has an LDSC-SEG result file -- this replaces the old
# ldsc/sumstats_urls.json lookup.
AG = Rn.set_index('trait')
AGM = pd.read_csv(f"{D}/bench_matrix.csv", index_col=0)
GRm = list(AGM.columns)

mrows = []
for e in panel:
    tr, exp, sid = e['trait'], e['expected'], e['study']
    respath = f"{P}/ldsc/results/{sid}.cell_type_results.txt"
    if not os.path.exists(respath):
        continue
    t = pd.read_csv(respath, sep='\t').sort_values('Coefficient_P_value').reset_index(drop=True)
    ld_rank = list(t.Name).index(exp) + 1
    ld_top = t.Name[0]
    ld_p = float(t.loc[t.Name == exp, 'Coefficient_P_value'].iloc[0])
    if tr in AGM.index:
        m = AGM.loc[tr]
        order = sorted(GRm, key=lambda gg: -m[gg])
        ag_rank = order.index(exp) + 1
        ag_top = order[0]
        ag_q = float(AG.loc[tr, 'q_exp'])
    else:
        ag_rank, ag_top, ag_q = np.nan, '--', np.nan
    mrows.append(dict(trait=tr, expected=exp, ag_rank=ag_rank, ag_top=ag_top, ag_q=ag_q,
                       ld_rank=ld_rank, ld_top=ld_top, ld_p=ld_p, agree=(ag_top == ld_top)))
C = pd.DataFrame(mrows)
C.to_csv(f"{D}/method_comparison.csv", index=False)

# ---- the "NMATCH-family" numbers (script lost; re-derived and verified
# against the values in the pre-existing numbers.tex -- see scripts/README.md)
NMATCH = len(C)
AGONE = int((C.ag_rank == 1).sum())
LDONE = int((C.ld_rank == 1).sum())
AGTHREE = int((C.ag_rank <= 3).sum())
LDTHREE = int((C.ld_rank <= 3).sum())
both_miss = C[(C.ag_rank > 1) & (C.ld_rank > 1)]
NBOTHMISS = len(both_miss)
NSAMEWRONG = int((both_miss.ag_top == both_miss.ld_top).sum())
NAGREE = int((C.ag_top == C.ld_top).sum())
PMATCH = sci(binomtest(AGONE, NMATCH, 1 / NT, alternative='greater').pvalue)
CHANCEMATCH = f"{NMATCH / NT:.1f}"
# NOVERLAP is dead code: it is defined in the historical numbers.tex but never
# used in main.tex (0 occurrences), and no combination of AG/LDSC rank or
# significance thresholds tried during reconstruction reproduces its stored
# value of 11 -- except the number of tissue groups itself. It is reproduced
# here as that (len(GR) == 11) for continuity; see scripts/README.md and the
# task report for the full account of what was tried.
NOVERLAP = len(GRm)

M['NMATCH'] = NMATCH
M['AGONE'] = AGONE
M['LDONE'] = LDONE
M['AGTHREE'] = AGTHREE
M['LDTHREE'] = LDTHREE
M['PMATCH'] = PMATCH
M['CHANCEMATCH'] = CHANCEMATCH
M['NBOTHMISS'] = NBOTHMISS
M['NSAMEWRONG'] = NSAMEWRONG
M['NAGREE'] = NAGREE
M['NOVERLAP'] = NOVERLAP
t2drep = C[C.trait == 'Type 2 diabetes (replication)']
if len(t2drep):
    r = t2drep.iloc[0]
    M['TWODLDP'] = f"{r.ld_p:.3f}"
    M['TWODAGRANK'] = int(r.ag_rank)
uc = C[C.trait.isin(['Ulcerative colitis (PheCode)', 'Ulcerative colitis'])]
if len(uc):
    r = uc.iloc[0]
    M['UCAGRANK'] = int(r.ag_rank)
    M['UCLDRANK'] = int(r.ld_rank)
    M['UCLDP'] = f"{r.ld_p:.3f}"
# per-trait matched-comparison macros for traits named in the prose (defined only when the trait is in the matched set)
for pref, tr in [('BMI', 'Body mass index'), ('COG', COGT), ('TOD', 'Type 1 diabetes'), ('UCX', 'Ulcerative colitis')]:
    sub = C[C.trait == tr]
    if len(sub):
        r = sub.iloc[0]
        _m, _e = f"{r.ld_p:.1e}".split("e")
        M[pref + 'LDRANK'] = int(r.ld_rank); M[pref + 'LDP'] = f"{r.ld_p:.3f}" if r.ld_p >= 0.001 else f"{_m}\\times 10^{{{int(_e)}}}"
        M[pref + 'LDTOP'] = NICE.get(r.ld_top, r.ld_top) if 'NICE' in globals() else r.ld_top
        M[pref + 'AGRANK'] = int(r.ag_rank)

open(f"{P}/numbers.tex", "w").write(
    "\n".join(f"\\newcommand{{\\{k}}}{{{v}}}" for k, v in M.items()) + "\n")

# ---- ldsc_table.tex: matched head-to-head table.
# Row order: stable sort by (LDSC rank, AlphaGenome rank) ascending, matching
# the historical table exactly (verified against data/v1_snapshot/ldsc_table.tex).
# Bold marks a rank-1 recovery in EITHER rank column independently -- this
# table has no trait-name-based bolding (unlike 31_tables.py / script 37's
# results_table.tex, which always bolds "Type 2 ..." rows regardless of rank;
# that convention is deliberately not applied here).
# Trait names are written RAW (as in bench_final.csv); the display-label
# relabelling for Type 1 diabetes / Ulcerative colitis is applied afterwards
# by scripts/37_review_analyses.py, exactly as it always was.
Cs = C.sort_values(['ld_rank', 'ag_rank'], kind='stable')


def b(v, ok):
    return f"\\textbf{{{v}}}" if ok else f"{v}"


lines = [
    "\\begin{table}[htbp]\\centering\\small",
    "\\caption{Matched head-to-head. Both methods were run on the same \\NMATCH{} traits, "
    "the same GWAS datasets and the same \\NTISSUE{} tissue definitions. \\textbf{Rank} is "
    "the position of the pre-specified expected tissue. Bold marks a rank-1 recovery.}\\label{tab:ldsc}",
    "\\begin{tabular}{llrlrl}\\toprule",
    "Trait & Expected & \\multicolumn{2}{c}{AlphaGenome} & \\multicolumn{2}{c}{LDSC-SEG} \\\\",
    "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}",
    " & & Rank & Top & Rank & Top \\\\ \\midrule",
]
for r in Cs.itertuples():
    lines.append(
        f"{esc(r.trait)} & {NICE[r.expected]} & {b(int(r.ag_rank), r.ag_rank == 1)} & "
        f"{NICE[r.ag_top]} & {b(int(r.ld_rank), r.ld_rank == 1)} & {NICE[r.ld_top]} \\\\")
lines += [
    "\\midrule",
    f"\\multicolumn{{2}}{{l}}{{\\emph{{rank-1 recovery}}}} & \\multicolumn{{2}}{{c}}{{{AGONE}/{NMATCH}}} "
    f"& \\multicolumn{{2}}{{c}}{{{LDONE}/{NMATCH}}} \\\\",
    f"\\multicolumn{{2}}{{l}}{{\\emph{{top-3 recovery}}}} & \\multicolumn{{2}}{{c}}{{{AGTHREE}/{NMATCH}}} "
    f"& \\multicolumn{{2}}{{c}}{{{LDTHREE}/{NMATCH}}} \\\\",
    "\\bottomrule\\end{tabular}\\end{table}",
]
open(f"{P}/ldsc_table.tex", "w").write("\n".join(lines))

print(f"n={n} traits scored | rank-1 {h1}/{n} | matched (LDSC) {NMATCH} | "
      f"AG rank-1 {AGONE}/{NMATCH} | LDSC rank-1 {LDONE}/{NMATCH}")
print("wrote: bench_summary.csv, bench_matrix.csv, bench_explain.csv, bench_final.csv, "
      "method_comparison.csv, ldsc_table.tex, numbers.tex")
