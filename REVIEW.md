# Internal review: "Tissue localisation by a sequence-to-function model tracks variant effect detectability, not tissue assay depth"

## Revision status (evening of 10 September 2026)

The author's decisions on the first review were: soften the title, replace the two underpowered datasets, run the brain head-to-head, and deposit the code and data on GitHub. Status:

- **Title** now reads "Per-variant tissue assignment by a sequence-to-function model is limited by effect detectability more than by assay coverage".
- **Panel revision.** Four datasets were replaced (all chosen as the largest GWAS of the trait with summary statistics both ingested by Open Targets and available from the GWAS Catalog): ulcerative colitis → Liu 2023 (GCST90446794, 23,252 cases); type 1 diabetes → Robertson 2021 (GCST90013445, 22,153 cases); body mass index → UK Biobank WGS 2025 (GCST90474606, 456,892 people); cognition → Savage 2018 intelligence (GCST006250, 269,867 people). The superseded datasets and their results are kept in a new Appendix C so nothing is hidden. Scoring used the same script logic as before (`scripts/38_score_panel_v2.py`), 516 variants, no failures.
- **Results on the revised panel** (`data/panel_v2.json`): expected tissue ranked first for 11 of 20 pairs (P = 2.7e-7), top three for 16, significantly enriched for 14. Both immune replacements recover immune tissue with significant enrichment (ulcerative colitis z = 0.50, P = 0.004; type 1 diabetes z = 0.82, P = 0.03), where the small PheCode datasets they replaced showed nothing. Both brain replacements still miss: BMI ranks brain 3rd (z = 0.29, P = 0.19), intelligence 3rd (z = 0.15, P = 0.87). The coverage correlations are now essentially zero (rho = -0.01 and +0.07; CI on the first -0.45 to +0.43). Visibility still does not predict rank across traits (rho = 0.31, P = 0.19); the within-trait per-variant effect remains (top-quartile odds ratio 1.45, P = 0.002; 1.78-fold excess versus the calibrated null in the top quartile).
- **Pipeline.** `scripts/39_rebuild.py <panel.json> .` now regenerates every benchmark file and `numbers.tex` from a panel definition, validated to reproduce the original outputs byte-for-byte on `data/panel_v1.json`; `scripts/37_review_analyses.py scores . <panel.json>` regenerates the review-stage numbers and tables. See `scripts/README.md`.
- **GitHub.** https://github.com/ronneldavis/alphagenome-tissue-benchmark (public). The current revision will be pushed once the LDSC runs are in.
- **LDSC-SEG for the four new datasets: done**, on the same GWAS and the same 11 tissue definitions (`ldsc/RUN_V2.md` has URLs, column choices, SNP counts and timings). Results: ulcerative colitis, immune first (P = 0.0017); type 1 diabetes, lung first and immune second (P = 0.016; bulk lung carries an immune expression signature); body mass index, brain first (P = 0.0016); intelligence, brain first (P = 2.4e-4). The matched comparison now covers 16 traits: LDSC-SEG ranks the expected tissue first for 10, AlphaGenome for 8; top three 15 versus 13; rank-1 and significant 10 versus 7. Both methods miss on 4 traits and pick the same wrong tissue on 2 (HDL cholesterol, atopic dermatitis).
- **What this means for the paper.** This was the test the thesis needed, and it passed: for the three traits where AlphaGenome misses and whose variants are least visible to it (BMI, intelligence, T2D), LDSC-SEG recovers the expected tissue from the same GWAS. The tissue signal is in the data; the per-variant predictions do not carry it. The abstract, introduction, Sections 2.5 and 2.6 and the Discussion now say this, and "Neither method dominates" became "LDSC-SEG is ahead", which the numbers require. The one case in AlphaGenome's favour remains coronary artery disease.
- **Abstract** rewritten and shortened to roughly 300 words. Cut further once you know the target journal's limit.

### Decisions still open
1. Target journal (sets the abstract limit and the format).
2. License for the repository (README says "to be decided").
3. Zenodo DOI for the deposit, to cite in Data availability, once you archive the GitHub repository.

The review of the first version follows unchanged.

## Review of the first version

Reviewed 10 September 2026 against the LaTeX source, the data in `data/`, the LDSC-SEG result files in `ldsc/results/`, the raw AlphaGenome score matrices (now copied to `scores/`), the GWAS Catalog and Open Targets APIs, and CrossRef.

## Verdict

The benchmark, the LDSC-SEG comparison and every number in the tables reproduce from the data. The bibliography is real: 24 of 25 entries were exactly right, one entry merged two papers, one author name was misspelled. The manuscript's problem was interpretation. Its central claim, that visibility predicts tissue recovery while assay depth does not, was overstated in two ways: the headline correlation (rho = 0.87) is largely a by-product of how the two quantities are computed, and the coverage test has too little power to support "not assay depth". Several factual statements about the inputs were also wrong (fine-mapping method, three phenotype definitions, summary-statistics availability), and one showcase example (ulcerative colitis) rests on a GWAS with 2,312 cases.

I have rewritten the manuscript so that every claim matches the data, added the analyses a referee would ask for, corrected the inputs' descriptions, and fixed the references. The revised PDF compiles cleanly. Six decisions remain for you before submission; they are listed at the end.

## What I changed in the manuscript

All edits are in `main.tex`, `refs.bib`, and the generated tables. Originals are kept as `*.orig` and `main.orig.pdf`. New numbers come from `scripts/37_review_analyses.py`, which writes `numbers_extra.tex` and regenerates `results_table.tex`, `traits_table.tex`, `tracks_table.tex`, the heatmap, and `data/pervariant_summary.csv`.

1. **Abstract rewritten** to state the shared-component caveat, the rank-based null result, the calibrated per-variant result, the significance-qualified LDSC comparison, and coronary artery disease (not ulcerative colitis) as the AlphaGenome-wins example.
2. **Section 2.4 (visibility) rewritten.** It now reports: inter-tissue correlation of predicted effects 0.80 (first component 75% of variance); visibility correlates with the top tissue's enrichment (rho = 0.88) as strongly as with the expected tissue's; visibility does not significantly predict rank across traits (rho = 0.29, P = 0.21); but at the variant level the expected tissue ranks first 1.81 times the calibrated null in the top visibility quartile versus 1.09 to 1.39 elsewhere, within-trait odds ratio 1.61 (95% CI 1.27 to 2.04, P = 8e-5).
3. **Section 2.3 (assay depth) softened** from "It does not" to "no evidence", with the 95% CI on rho (-0.63 to +0.21), the rank-based versions, and a note that the islet example is non-significant and not replicated.
4. **Section 2.1 corrected**: fine-mapping is SuSiE for 3 studies, PICS on summary statistics for 10, and PICS around reported lead SNPs for 7 (so those 7 score lead SNPs, not fine-mapped variants); the background is one shared set near T2D and LDL loci, not matched on frequency, LD or annotation.
5. **Section 2.2**: results table gains a bootstrap-support column (4 of the 10 rank-1 recoveries have support below 0.5); the unique-trait and calibrated-chance versions of the binomial test are given; ulcerative colitis is explicitly not counted as a success.
6. **Section 2.5 (LDSC)**: "Neither method wins" became "Neither method dominates", with the significance-qualified count (LDSC-SEG 7, AlphaGenome 5 of 14); coronary artery disease replaces ulcerative colitis as the reverse example; atopic dermatitis is named as the second same-wrong-tissue case (LDSC immune P = 2.4e-5); the summary-statistics criterion is stated correctly.
7. **Section 2.6 (T2D)**: the Open Targets caveat now describes what actually happened (PICS around reported index variants, no summary statistics ingested); the "same variants" wording is fixed (LDSC-SEG uses the whole GWAS); the LDSC P is flagged as nominal with bulk pancreas as the islet proxy.
8. **Discussion**: "roughly a third" became 8 of 20; polygenicity argument now admits platelet count as a counterexample; two alternative readings of low visibility are added (non-causal scored variant; marginal GWAS signal), with the supporting correlations; the false "no public summary statistics" statement is replaced.
9. **Limitations**: three new items (heterogeneous fine-mapping, unmatched background, low power for the coverage test); the proxy-phenotype item now names the actual phenotypes and case counts; item 9 no longer cites atopic dermatitis as a non-significant LDSC case (it is highly significant for immune); "revised two labels" became "contradicts two labels".
10. **Methods**: fine-mapping regimes, coding-variant handling (signals dropped, not replaced), background construction, bootstrap, calibrated null, per-variant models and inter-tissue correlation are described.
11. **Tables**: trait table gains GWAS sample size, case counts, and fine-mapping regime, and its variant-count column is no longer labelled N; track census gains a five-mark histone column (the count actually used in the correlations, 45 for brain, not 63); three trait labels corrected everywhere (tables, heatmap, LDSC table).
12. **Bibliography**: `sasse2023personal` split into `sasse2023benchmark` (Sasse et al., Nat Genet 55:2060, doi 10.1038/s41588-023-01524-6) and `huang2023personal` (Huang et al., Nat Genet 55:2056, doi 10.1038/s41588-023-01574-w), both now cited; Bonàs-Guarch spelling fixed; Maurano et al. 2012 added for the non-coding claim.

## Findings in detail

### Major (would block acceptance)

**1. The visibility correlation is largely mechanical.** Visibility is the mean standardised effect over all 294 tissue-assigned histone tracks; the recovery measure it was correlated with is the expected tissue's mean standardised effect, which is a subset of the same tracks. Predicted histone effects are dominated by a shared component: mean pairwise inter-tissue Spearman correlation across the 2,614 scored variants is 0.80 (range 0.68 to 0.91) and the first principal component carries 75% of the variance. So visibility correlates with the top tissue's enrichment, whichever tissue that is, at rho = 0.88. The scale-free test (does visibility predict the expected tissue's rank?) gives rho = 0.29, P = 0.21; bootstrap support rho = 0.31, P = 0.18; median visibility for rank-1 versus other pairs 0.69 versus 0.31, P = 0.12. Six of the ten most visible pairs rank the expected tissue first, and so do four of the ten least visible. The paper's sentence "Where it does [predict anything], the tissue it names is usually right" was not supported. Fixed in the text; the per-variant analysis (below) is the salvage.

**2. Per-variant test, properly calibrated.** The per-variant null is not uniform: after per-tissue standardisation, immune ranks first for 30.5% of background variants and liver for 5.0%. Against a null calibrated per expected tissue, the expected tissue ranks first for 16.0% of variants versus 11.8% expected (1.36-fold, P = 1e-10). The excess concentrates in the top visibility quartile (1.81-fold, P = 5e-12) versus 1.09 to 1.39 in the other three. A logistic regression with trait fixed effects gives an odds ratio of 1.61 for a trait's own top-quartile variants (P = 8e-5) and 1.12 per SD of visibility (P = 0.02). This is real support for the detectability account at the variant level. It is modest, and the paper now says so.

**3. The coverage test cannot support "not assay depth".** Track count takes eight distinct values across the 20 pairs. The 95% CI on rho(track count, enrichment) is -0.63 to +0.21. Rank-based versions are rho = +0.09 (P = 0.69) and +0.15 (P = 0.53). The islet extreme is a non-significant rank-1 (P = 0.07, bootstrap support 0.44) that reverses in the replication GWAS. The text now says the data exclude a strong coverage effect, not a moderate one. The title still asserts "not tissue assay depth"; see decisions.

**4. Fine-mapping was misdescribed.** The paper said Open Targets "applies SuSiE-style fine mapping to GWAS summary statistics". The Platform's own records (`data/finemapping_provenance.json`) show: SuSiE for 3 studies (FinnGen Crohn's, CAD 98%, Xue T2D 93%), PICS on summary statistics for 10, PICS around reported top hits with no summary statistics for 7 (Mahajan T2D, triglycerides, asthma, platelet count, BMI, cognition, lung function). For the last group the top-PIP variant is the reported lead SNP. This also invalidated the T2D caveat paragraph ("Open Targets re-derives credible sets from summary statistics") and changes the interpretation of the 15-variant overlap between the two T2D sets (different GWAS and different methods).

**5. Three phenotypes are not what the paper calls them.**
- "Type 1 diabetes (PheCode)" is GCST90479877, PheCode 250.13 "Type 1 diabetes with ophthalmic manifestations", 3,569 cases. Its 72 credible sets are 94% singletons with PIP 1.0 at median -log10 P of 9. Its rank-7 failure says nothing about the model.
- "Ulcerative colitis (PheCode)" is GCST90480318, PheCode 555.21, 2,312 cases; 149 credible sets, 98% singletons with PIP 1.0, median -log10 P of 9. A GWAS this small does not yield 149 real loci. The abstract used this trait as AlphaGenome's showcase win over LDSC-SEG (immune ranked first); its immune z is -0.01. Removed as a success everywhere.
- "Crohn's disease (medication proxy)" is FinnGen RX_CROHN_1STLINE, "first line medication for Crohn's disease", 106,320 cases out of 500,348. One in five Finns is not a Crohn's case; this is a broad medication-purchase endpoint.
- "Cognitive ability / education" is GCST008595, Lam et al. 2019, a pleiotropic meta-analysis of cognition, educational attainment and schizophrenia, retrieved under the schizophrenia ontology term (`data/trait_panel.json`). Brain remains the right expected tissue, but the label was misleading.
All four are relabelled in tables, figure and text, with case counts in the trait table.

**6. "No public genome-wide summary statistics" for BMI and cognition is false.** The specific GWAS used (Kichaev 2019, Lam 2019) are not distributed through the GWAS Catalog FTP, which was the pipeline's actual criterion (`ldsc/scripts/90_pipeline.py`). BMI summary statistics from GIANT (Locke 2015, Yengo 2018) and cognition or educational attainment statistics (Savage 2018, Lee 2018) are public. A referee will ask why the strongest test of the paper's thesis was not run. The LDSC environment from the earlier session is intact at the path in `ldsc/scripts/tissues.ldcts`, so this is hours, not days.

**7. "Neither method wins" understated LDSC-SEG.** Requiring rank 1 and significance for the expected tissue: LDSC-SEG 7 of 14, AlphaGenome 5 of 14 (ulcerative colitis q = 0.99 and atrial fibrillation q = 0.17 are AlphaGenome rank-1 calls without enrichment). LDSC-SEG also leads on top-three, 12 versus 10. Text corrected.

**8. Limitation 9 was wrong about atopic dermatitis.** LDSC-SEG gives immune P = 2.4e-5 for atopic dermatitis, the strongest tissue call in the whole LDSC run, with skin at rank 7 (P = 0.59). It is not a "non-significant ordering" case; it is the second trait on which both methods confidently reject the pre-specified label. This strengthens the "labels are wrong" argument and the text now uses it. Coronary artery disease is the genuine no-tissue-significant LDSC case (best P = 0.053).

**9. Background is not "matched for genomic context".** `scripts/11_background.py` draws two positions 50 to 400 kb from each non-coding T2D and LDL signal, keeps 320 at random, and converts to transition SNVs. There is no matching on allele frequency, LD, gene density or annotation, and the same 318 variants serve every trait. Enrichment over this background partly measures that GWAS variants are regulatory at all; rank is the tissue-specific readout. The abstract and methods now say this.

### Moderate

- **20 pairs, 19 traits.** T2D appears twice. The binomial test treats pairs as independent. Counting T2D once: P = 1.6e-6 (Mahajan kept) or 1.7e-5 (Xue kept). Both versions are now given.
- **Calibrated chance.** Resampling background sets of 150 variants gives an expected 1.7 rank-1 recoveries rather than 1.8, because tissue tails differ. Now stated.
- **Track counts inconsistent.** The text quoted brain's 63 histone tracks; the correlation used the five-mark count, 45. Both are now in the track table.
- **"Roughly a third" fall below threshold** was 8 of 20 (40%). Fixed.
- **Polygenicity argument** ignored platelet count, which is highly polygenic and recovered strongly. Now acknowledged, with two alternative explanations for low visibility that the design cannot exclude: the scored variant is not causal (visibility rises with PIP, rho = 0.11 per variant, P = 6e-9; with GWAS significance rho = 0.09, P = 1e-6; trait-level median -log10 P versus expected-tissue enrichment rho = 0.47, P = 0.04), and the input GWAS is underpowered.
- **Coding-variant handling.** When the top-PIP variant is coding the signal is dropped, not replaced by the next non-coding member (`scripts/22_run_all.py`). Methods now say so.
- **LDSC T2D P = 0.032** is nominal across 11 tissues and uses bulk pancreas. Now flagged at the point of use.
- **Rank-1 fragility.** Bootstrap support: T2D 0.44, Crohn's proxy 0.48, ulcerative colitis 0.49, platelet count 0.44, heart rate 0.63; the liver traits, atrial fibrillation and CAD are above 0.75. Heart rate is instructive: heart ranks first at the trait level but for individual variants only 8.7% versus 8.5% null, so the call rests on a few large-effect loci.

### Minor

- Uncited bibliography entries (`miguelescalada20193d`, `grant2006tcf7l2`, `prokopenko2009mtnr1b`, `dupuis2010fasting`, `aragam2022cad`, `fulco2019abc`) are harmless with natbib but could be pruned. Optional issue numbers for several entries are listed in the reference table below.
- The sentence "an error we made and corrected during this work" (Section 2.1) is unusual in a paper. It is honest; keep or cut as you prefer.
- `\date{\today}` will change at every compile.
- Data availability says "available from the author on request". Most journals now require deposition (Zenodo or GitHub for scripts and score matrices). The `scores/` folder (13 MB) plus `data/` and `ldsc/results/` is everything needed.
- The acknowledgement of Claude Code is appropriate; check the target journal's policy wording on AI assistance.

## Reference check (CrossRef, all 25 entries)

| Key | Status | Note |
|---|---|---|
| avsec2026alphagenome | correct | Nature 649(8099):1206-1218, 2026; DOI resolves; preprint 10.1101/2025.06.25.661532 |
| avsec2021enformer | correct | issue 10 optional |
| linder2025borzoi | correct | pages 949-961 confirmed, issue 4 |
| sasse2023personal | **wrong, fixed** | title and DOI were Huang et al.; authors and pages were Sasse et al. Split into two entries |
| karollus2023promoters | correct | |
| finucane2018seg | correct | |
| finucane2015partitioning | correct | |
| chiou2021islet | correct | published title uses an en dash "cell type- and state-specific" |
| gaulton2010islet | correct | |
| pasquali2014enhancer | correct | |
| miguelescalada20193d | **author fixed** | Bonàs-Guarch, not Bonhs-Guarch (uncited) |
| greenwald2019islet | correct | |
| udler2018clusters | correct | |
| mahajan2018finemap | correct | |
| xue2018t2d | correct | |
| grant2006tcf7l2 | correct | uncited |
| prokopenko2009mtnr1b | correct | uncited |
| dupuis2010fasting | correct | uncited |
| buniello2025opentargets | correct | |
| gtex2020 | correct | |
| encode2020 | correct | |
| wang2020susie | correct | |
| aragam2022cad | correct | uncited |
| lotta2017ir | correct | |
| fulco2019abc | correct | uncited |
| maurano2012regulatory | **added** | Science 337(6099):1190-1195 |

GWAS accessions: all 19 GCST accessions and the FinnGen endpoint resolve, and first authors and years in the trait table match the GWAS Catalog.

## Numerical checks

Every value in `numbers.tex`, `results_table.tex`, `ldsc_table.tex` and `tracks_table.tex` was recomputed from `data/bench_final.csv`, `data/method_comparison.csv`, `ldsc/results/*.txt` and the track CSVs. All matched: 10/20 rank-1 (P = 3.0e-6), 15/20 top-3 (P = 1.2e-5), 12 BH-significant, rho values, MWU P, 7/7 and 10/12 in the matched comparison, 5 both-miss and 2 same-wrong, LDSC ranks and P-values, histone and accessibility track counts. The only errors were in prose, not in the tables.

## Decisions for the author

1. **Title.** "tracks variant effect detectability, not tissue assay depth" now overstates both halves. Candidates: "Tissue localisation of GWAS variants by AlphaGenome: recovery above chance, a dominant tissue-shared effect component, and parity with LDSC-SEG", or "Per-variant tissue assignment by AlphaGenome is limited by effect detectability more than by assay coverage". Your call; I left the title unchanged.
2. **Underpowered traits.** Consider replacing the two MVP PheCode traits (ulcerative colitis 2,312 cases; T1D-ophthalmic 3,569 cases) with de Lange 2017 IBD and Chiou 2021 or Robertson 2021 T1D, which Open Targets fine-maps with SuSiE. This needs an AlphaGenome API rerun (about 300 variants); no key was available in this environment.
3. **Brain traits head-to-head.** Run both methods on a BMI and a cognition GWAS that have public summary statistics. This is the single most persuasive addition possible and the one a referee will insist on.
4. **The new analyses.** Section 2.4, the bootstrap column, the calibrated null and the new Methods text are mine. Please read them as critically as a referee would; the code is `scripts/37_review_analyses.py`, and every number in the text is a macro from `numbers_extra.tex`.
5. **Deposition.** Put `scores/`, `data/`, `ldsc/results/`, `scripts/` and `ldsc/scripts/` on Zenodo or GitHub and cite the DOI in Data availability.
6. **Housekeeping.** The abstract is now about 520 words, which is over most journals' limit (typically 150 to 300); it needs cutting once you have settled the framing. Also: uncited bib entries; the "error we made" sentence; the `\today` date; `\CODEURL` is defined but unused.

## What I did not do

No new AlphaGenome scoring (no API key here) and no new LDSC-SEG runs. All added statistics reuse the existing score matrices.

## Reproduction

```bash
uv run --with numpy --with pandas --with scipy --with statsmodels --with matplotlib \
    python scripts/37_review_analyses.py scores .
tectonic main.tex
```

`scripts/36_numbers2.py` still writes `numbers.tex`; the new macros live separately in `numbers_extra.tex` so the two do not overwrite each other.
