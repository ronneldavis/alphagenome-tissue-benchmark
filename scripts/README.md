# Scripts

How to rerun the benchmark and regenerate the paper's tables and numbers, and
what every script in this directory (and `ldsc/scripts/`) did on the way to
the current data.

## Reproduction order

```bash
uv run --with numpy --with pandas --with scipy --with statsmodels --with matplotlib \
    python scripts/39_rebuild.py data/panel_v1.json .

uv run --with numpy --with pandas --with scipy --with statsmodels --with matplotlib \
    python scripts/37_review_analyses.py scores . data/panel_v1.json

tectonic main.tex
```

1. **`39_rebuild.py <panel.json> <P>`** runs the core benchmark from the
   score files in `scores/` and the data CSVs in `data/`: per-tissue
   background standardisation, rank and significance of the expected tissue,
   visibility, track-coverage counts, and the matched AlphaGenome-vs-LDSC-SEG
   comparison. Writes `data/bench_summary.csv`, `data/bench_matrix.csv`,
   `data/bench_explain.csv`, `data/bench_final.csv`,
   `data/method_comparison.csv`, `ldsc_table.tex` and `numbers.tex`.
2. **`37_review_analyses.py <scores> <P> <panel.json>`** adds the review-stage
   analyses (bootstrap support, the calibrated per-variant null, the
   shared-component check, fine-mapping regime) on top of `39_rebuild.py`'s
   output. It also relabels a few trait names to their display form in
   `ldsc_table.tex` (in place) and regenerates `results_table.tex`,
   `traits_table.tex`, `tracks_table.tex`, the heatmap and
   `data/pervariant_summary.csv`. Run it *after* `39_rebuild.py`, never
   before -- it reads `data/bench_final.csv` and the existing `ldsc_table.tex`.
3. **`tectonic main.tex`** compiles the paper from those `.tex` files. Not run
   by either script above.

Both scripts take the panel file as an explicit argument now; nothing in
either one is hard-coded to the current 20 trait-GWAS pairs. `data/panel_v1.json`
is the panel behind every number currently in the paper. `data/panel_v2.json`
swaps in four better-powered replacement studies (see "Panel versions" below)
and is not yet used anywhere -- its score files exist under `scores/` but
`39_rebuild.py` has only been run, and validated, against `panel_v1.json`.

### Two things that will bite you if you skip this paragraph

- **Panel row order is part of the reproducibility contract for
  `37_review_analyses.py`.** Its bootstrap (2,000 resamples per trait) draws
  from one `numpy` random generator shared across every trait, in the
  panel file's row order. Reorder an existing panel file's rows and every
  bootstrap number from the moved row onward changes -- not a bug, just a
  different, equally valid draw from the same seed. This is why
  `data/panel_v1.json` lists "Type 2 diabetes (replication)" *last*: that
  matches the row order the old hard-coded `FILES` dict in
  `37_review_analyses.py` used, and reordering it here (to match
  `30_bench.py`'s old `LABEL` dict order instead) silently changed the
  bootstrap column until it was moved back. `39_rebuild.py` has no such
  sensitivity -- its outputs are aggregate counts and sorts, not a shared
  random stream -- so this only matters for the panel file's row order,
  not for anything else about it.
- **`data/accessions.json` is keyed by trait display name, not by panel
  version.** `panel_v2.json` reuses the name "Body mass index" for a new
  study (GCST90474606) while `panel_v1.json`'s "Body mass index" is still
  GCST007039. Both cannot occupy the same JSON key, and overwriting it would
  silently break `traits_table.tex` the next time `panel_v1.json` is
  rebuilt. `accessions.json`'s "Body mass index" entry was therefore left at
  GCST007039; if `panel_v2.json` is ever actually run through
  `37_review_analyses.py`, that one entry needs a manual decision (update it,
  knowing it then reads wrong for `panel_v1.json`; or give the v2 trait a
  distinguishing label first).

### The lost "NMATCH-family" macros

`numbers.tex` carries a block of macros (`NMATCH`, `AGONE`, `LDONE`,
`AGTHREE`, `LDTHREE`, `PMATCH`, `CHANCEMATCH`, `NBOTHMISS`, `NSAMEWRONG`,
`TWODLDP`, `TWODAGRANK`, `UCAGRANK`, `UCLDRANK`, `NAGREE`, `NOVERLAP`) that
summarise the matched AlphaGenome-vs-LDSC-SEG comparison. No script producing
them survived in `scripts/`; `39_rebuild.py` re-derives all fifteen from
`data/method_comparison.csv`, and every one except `NOVERLAP` was checked
against the historical `numbers.tex` value and matches exactly. `NOVERLAP`
is dead: it is never referenced in `main.tex`, and no rank- or
significance-based formula tried during reconstruction reproduces its stored
value of 11 -- except the plain number of tissue groups, which is what
`39_rebuild.py` now computes for it. It does not affect anything in the
compiled paper.

## Panel versions

- **`data/panel_v1.json`** -- the 20 trait-GWAS pairs behind the current
  paper. One list entry per pair: `trait` (display name, matches
  `data/bench_final.csv`), `label` (table/figure display label, usually
  equal to `trait`), `study` (GWAS Catalog / FinnGen accession), `expected`
  (pre-specified tissue group), `source` (first-author citation shown in
  tables), `file` (score-file prefix under `scores/`).
- **`data/panel_v2.json`** -- the same panel with four replacements chosen to
  address REVIEW.md's decisions 2 and 3: the two underpowered MVP PheCode
  traits (ulcerative colitis, 2,312 cases; type 1 diabetes, 3,569 cases) swap
  for larger, SuSiE- or PICS-on-summary-statistics-fine-mapped GWAS (Liu 2023
  ulcerative colitis; Robertson 2021 type 1 diabetes), and BMI and cognition
  swap for GWAS with public summary statistics (UKB-WGS 2025 BMI; Savage 2018
  intelligence), enabling the LDSC-SEG head-to-head that the original panel
  could not run for those two. All four new studies have `hasSumstats: true`
  in Open Targets, confirmed when their entries were added to
  `data/study_meta.json` and `data/finemapping_provenance.json`.

## What the original scripts did

Early exploration, superseded by the numbered scripts that follow them:

- `01_tracks.py`, `02_tissues.py` -- listed AlphaGenome's available data
  types and tissue names, looking for pancreas/islet coverage.
- `03_panel.py`, `09_panel_def.py` -- early 5-tissue-group drafts, superseded
  by `20_panel.py`'s final 11-group version.
- `04_fetch_variants.py`, `05_build_variants.py` -- fetched and cleaned
  fine-mapped variants for one diabetes study; the single-trait prototype of
  `22_run_all.py`.
- `06_probe_score.py`, `07_probe_one.py` -- one-off checks of AlphaGenome's
  scoring API and output shape.
- `08_score.py` -- scored one trait's variants; superseded by
  `22_run_all.py`.
- `10_analyse.py` -- first tissue-enrichment comparison (diabetes vs.
  cholesterol, no background normalisation); superseded by `12_analyse2.py`.
- `12_analyse2.py` -- re-ran that comparison against the random background
  set built by `11_background.py`, the direct ancestor of `30_bench.py`'s
  method.
- `13_followup.py` -- looked at which individual genes/variants drove the
  liver and islet signals (including TCF7L2).
- `14_histone_only.py` -- reran the enrichment check restricted to histone
  tracks, to check it wasn't an artefact of which tracks were included.
- `15_loci.py`, `16_export.py`, `17_gen.py` -- per-variant tissue-specificity
  flags, a JSON export and an HTML report for the early 2-trait pilot;
  superseded by `40_artifact.py` / `41_write_art.py`.

The panel and scoring run behind the current data:

- `20_panel.py` -- defined the final 11 tissue groups and tagged every
  AlphaGenome track with one.
- `21_traits.py` -- searched Open Targets for the best-powered GWAS per
  trait, producing `data/trait_panel.json`.
- `22_run_all.py` -- fetched fine-mapped variants and scored them with
  AlphaGenome for every trait in the panel; this is what produced the
  `scores/T_<study>_*` files. Drops a signal instead of substituting the
  next credible-set member when its top-PIP variant is coding.
- `11_background.py` -- built the 318-variant shared background set (near
  T2D/LDL loci, not matched on frequency/LD/annotation) that every trait is
  compared against.
- `23_studymeta.py` -- fetched publication metadata (first author, date,
  journal, PMID, N) per study from Open Targets into `data/study_meta.json`.
  The same query pattern, extended with fine-mapping method and confidence,
  built `data/finemapping_provenance.json`.
- `38_score_panel_v2.py` -- the same fetch-and-score logic as
  `22_run_all.py`, restricted to the four `panel_v2.json` replacement
  studies, with a `SMOKE=1` single-study test mode.

The benchmark itself (see `39_rebuild.py`'s docstring for the exact formula
lineage):

- `30_bench.py` -- the core benchmark: per-tissue background standardisation
  over the five histone marks, rank of the expected tissue, one-sided
  Mann-Whitney P. Originally hard-coded the trait/study/tissue/source/file
  mapping now in the panel file.
- `33_t2d_deep.py` -- added the BH-corrected q-value (`q_exp`) and the
  Mahajan-vs-Xue T2D variant-overlap numbers; the q-value step has no
  surviving standalone script reference but is folded into `39_rebuild.py`.
- `34_explain.py` -- added the track-count columns (`n_tracks_exp`,
  `acc_tracks_exp`) reused by every later script; its own `visibility`
  column (`bench_matrix.max(axis=1)`) was mechanical (REVIEW.md finding 1)
  and is not reproduced anywhere in `39_rebuild.py`.
- `35_visibility.py` -- the real `visibility`: mean standardised effect
  pooled over every tissue-assigned five-mark track, independent of which
  tissue is expected. This is the version `39_rebuild.py` uses everywhere.
- `32_numbers.py` -- first pass at the `numbers.tex` macros; superseded by
  `36_numbers2.py`.
- `36_numbers2.py` -- the macros actually in the paper's `numbers.tex`
  (hit counts, binomial P, correlations, the two T2D datasets' variant
  overlap).
- `31_tables.py` -- the original `results_table.tex` / `traits_table.tex` /
  `tracks_table.tex` / heatmap, superseded for the current paper by the
  versions `37_review_analyses.py` writes (bootstrap column, GWAS size and
  case counts, five-mark histone column).
- `37_review_analyses.py` -- the review-stage additions: bootstrap support,
  the calibrated per-variant null and visibility quartiles, the inter-tissue
  shared-component check, fine-mapping regime per study, and the relabelled
  tables/heatmap/`ldsc_table.tex`. Writes `numbers_extra.tex` separately so
  it never collides with `36_numbers2.py`'s `numbers.tex`.
- `40_artifact.py`, `41_write_art.py` -- built chart images and an HTML
  write-up of the (pre-review) 20-trait findings; not part of the LaTeX
  build.

LDSC-SEG comparison (`ldsc/scripts/`, plus `scripts/50_ldsc_compare.py`):

- `60_genesets.py` -- built tissue-specific gene sets from GTEx v8 expression.
- `61_make_annot.py` -- annotated SNPs by proximity to each tissue's genes,
  per chromosome.
- `63_split.py` -- split per-SNP annotations into per-tissue files and wrote
  `tissues.ldcts`.
- `90_pipeline.py` -- the full run: downloaded and munged GWAS summary
  statistics, then ran LDSC-SEG per trait, producing `ldsc/results/*.txt`.
- `50_ldsc_compare.py` -- compared AlphaGenome's top tissue against
  *published* LDSC-SEG findings from Finucane et al. 2018 (not a fresh run);
  a sensitivity check, superseded for the paper's own comparison by
  `ldsc/scripts/95_compare.py`.
- `ldsc/scripts/95_compare.py` -- the matched head-to-head: same traits, same
  GWAS, same 11 tissue definitions, both methods. Originally read the
  trait/accession/tissue mapping from `ldsc/sumstats_urls.json`;
  `39_rebuild.py` reproduces its exact logic but takes that mapping from the
  panel file instead, and treats a trait as "matched" iff
  `ldsc/results/<study>.cell_type_results.txt` exists.
