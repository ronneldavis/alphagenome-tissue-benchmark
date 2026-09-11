# AlphaGenome panel v2 scoring — validation report

Job: `scripts/38_score_panel_v2.py scores`
Log: `scores/score_v2.log`

- Start: 2026-09-10 20:23:48 EDT (Infisical secret injection line)
- End: 2026-09-10 20:32:11 EDT (log shows `ALL DONE`, process exited)
- Total wall time: ~8m 23s
- Failures: none (`n_failed=0` for all four studies)

## Per-study validation

| Study | GCST | n scored | median pip | median negp | share cs_size==1 | NaN rows in hist |
|---|---|---|---|---|---|---|
| Ulcerative colitis (Liu 2023) | GCST90446794 | 150 | 0.9995 | 13.0 | 0.660 | 0 |
| Type 1 diabetes (Robertson 2021) | GCST90013445 | 66 | 0.2631 | 13.5 | 0.076 | 0 |
| Body mass index (UKB-WGS 2025) | GCST90474606 | 150 | 0.2490 | 16.0 | 0.087 | 0 |
| Intelligence (Savage 2018) | GCST006250 | 150 | 0.2531 | 10.0 | 0.073 | 0 |

All four studies passed schema checks:
- `hist` array shape is `(n, 1116)` for every study.
- `rna` array is 2-D with the same row count `n` as `hist` (371 columns each).
- `_meta.csv` has exactly `n` rows with columns `variant_id, rsid, pip, l2g_gene, cs_size`.
- `_variants.csv` has the required columns (`variant_id, rsid, chrom, pos, ref, alt, pip, cs_size, csq, l2g_gene, negp`) and row counts `>= n` (candidate variants before capping/filtering).
- All `n` scored `variant_id`s from `_meta.csv` were found in `_variants.csv` (used to pull `negp` for the median).

No anomalies found.
