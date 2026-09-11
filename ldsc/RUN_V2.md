# LDSC-SEG run v2 — 4 additional GWAS traits

This run adds four traits to the LDSC-SEG tissue-enrichment analysis, using the
patched LDSC fork and reference files already in place for the rest of the
pipeline (`baselineLD`, `1000G_Phase3_weights_hm3_no_MHC`, `w_hm3.snplist`,
`tissues.ldcts`, 11 tissue bins: ADIPOSE, ARTERY, BRAIN, HEART, IMMUNE,
INTESTINE, ISLET, LIVER, LUNG, MUSCLE, SKIN).

Scripts: `scripts/90b_pipeline_v2.py` (original driver — download, munge,
`--h2-cts`, for all four), `scripts/90c_finish.py` (finished the three
already-downloaded traits without re-downloading; also supports a one-trait
rerun mode used to apply a column fix), `scripts/90d_bmi.py` (finished the
fourth trait, BMI, once its download was repaired).

All four `--h2-cts` outputs have the expected 11 rows (one per tissue) with
columns `Name, Coefficient, Coefficient_std_error, Coefficient_P_value`.

## Segmented download note (body mass index)

The BMI file is 3.96 GB. Its single-connection resumable download (`curl -C -`)
slowed to about 190 KB/s and was going to take hours, so the run was stuck. The
job was killed (`SIGTERM` to the pipeline process and its `curl` child) without
deleting the partial file (3,671,203,840 of 3,956,491,232 bytes, 92.8% done).

The remaining 285,287,392 bytes were fetched as 8 parallel `curl -r START-END`
byte-range requests (~35.66 MB each, contiguous, starting right after the
partial file's last byte), each with its own retry loop. All 8 segments
completed on the first attempt — no retries were needed. They were appended to
the partial file in order, the final size was checked byte-for-byte against the
server's `Content-Length` (3,956,491,232 — exact match), and `gzip -t` passed.
Wall-clock time for the segmented portion was about 2–3 minutes, versus an
estimated multi-hour finish at the original single-connection rate.

## Per-study results

### Ulcerative colitis — GCST90446794 (expected tissue: IMMUNE)

- URL: `https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90446001-GCST90447000/GCST90446794/harmonised/GCST90446794.h.tsv.gz`
- File size: 822,519,117 bytes (822.5 MB) — verified against remote `Content-Length` and `gzip -t` before use
- Columns used: `--snp rsid --a1 effect_allele --a2 other_allele --p p_value --signed-sumstats beta,0 --N 375508`
- Munge: 1,190,321 SNPs, mean χ² = 1.467 (43.1 s)
- `--h2-cts` run time: 20 m 24.5 s (total incl. munge: ≈ 21 m 8 s)
- 11 tissues by `Coefficient_P_value` (ascending):

  | Rank | Tissue | P-value |
  |---|---|---|
  | 1 | IMMUNE | 0.0017 |
  | 2 | LUNG | 0.015 |
  | 3 | INTESTINE | 0.033 |
  | 4 | ARTERY | 0.43 |
  | 5 | ADIPOSE | 0.49 |
  | 6 | ISLET | 0.65 |
  | 7 | MUSCLE | 0.79 |
  | 8 | LIVER | 0.86 |
  | 9 | BRAIN | 0.92 |
  | 10 | HEART | 0.94 |
  | 11 | SKIN | 0.99 |

- **Expected tissue (IMMUNE) rank: 1 of 11.**

### Type 1 diabetes — GCST90013445 (expected tissue: IMMUNE)

- URL: `https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90013001-GCST90014000/GCST90013445/harmonised/GCST90013445.h.tsv.gz`
- File size: 12,886,325 bytes (12.9 MB) — verified against remote `Content-Length` and `gzip -t` before use
- Columns used: `--snp rsid --a1 effect_allele --a2 other_allele --p p_value --signed-sumstats beta,0 --N 59527`
- Munge: 1,190,321 SNPs, mean χ² = 6.473 (10.4 s)
- `--h2-cts` run time: 13 m 53.0 s (total incl. munge: ≈ 14 m 3 s)
- 11 tissues by `Coefficient_P_value` (ascending):

  | Rank | Tissue | P-value |
  |---|---|---|
  | 1 | LUNG | 0.012 |
  | 2 | IMMUNE | 0.016 |
  | 3 | ISLET | 0.14 |
  | 4 | INTESTINE | 0.19 |
  | 5 | LIVER | 0.27 |
  | 6 | BRAIN | 0.33 |
  | 7 | MUSCLE | 0.72 |
  | 8 | ARTERY | 0.77 |
  | 9 | HEART | 0.87 |
  | 10 | SKIN | 0.89 |
  | 11 | ADIPOSE | 0.92 |

- **Expected tissue (IMMUNE) rank: 2 of 11** (P = 0.016, just behind LUNG at P = 0.012).

### Body mass index — GCST90474606 (expected tissue: BRAIN)

- URL: `https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90474001-GCST90475000/GCST90474606/harmonised/GCST90474606.h.tsv.gz`
- File size: 3,956,491,232 bytes (3.96 GB) — finished via segmented download (see note above); final size verified byte-exact against remote `Content-Length` and `gzip -t` passed
- Columns used: `--snp rsid --a1 effect_allele --a2 other_allele --p p_value --signed-sumstats beta,0 --N-col n`
- Munge: 1,190,321 SNPs, mean χ² = 3.309 (4 m 48.7 s — this is the largest raw input of the four)
- `--h2-cts` run time: 13 m 42.0 s (total incl. munge: ≈ 18 m 31 s; excludes download time)
- 11 tissues by `Coefficient_P_value` (ascending):

  | Rank | Tissue | P-value |
  |---|---|---|
  | 1 | BRAIN | 0.0016 |
  | 2 | ISLET | 0.61 |
  | 3 | ARTERY | 0.90 |
  | 4 | INTESTINE | 0.95 |
  | 5 | MUSCLE | 0.96 |
  | 6 | SKIN | 0.97 |
  | 7 | LIVER | 0.97 |
  | 8 | ADIPOSE | 0.99 |
  | 9 | LUNG | 1.00 |
  | 10 | IMMUNE | 1.00 |
  | 11 | HEART | 1.00 |

- **Expected tissue (BRAIN) rank: 1 of 11.**

### Intelligence — GCST006250 (expected tissue: BRAIN)

- URL: `https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006250/harmonised/29942086-GCST006250-EFO_0004337.h.tsv.gz`
- File size: 570,198,602 bytes (570.2 MB) — verified against remote `Content-Length` and `gzip -t` before use
- Columns used: `--snp hm_rsid --a1 hm_effect_allele --a2 hm_other_allele --p p_value --signed-sumstats hm_beta,0 --N-col n_analyzed --ignore effect_allele,other_allele`
- **Column-name fix:** the first munge attempt failed (`AttributeError: 'DataFrame' object has no attribute 'str'`, 0 SNPs). Cause: this harmonised file carries both the `hm_`-prefixed columns we asked for and the original unprefixed `effect_allele`/`other_allele` columns. `munge_sumstats.py`'s built-in column-name dictionary auto-matched the unprefixed pair to A1/A2 as well, on top of our explicit `--a1`/`--a2`, producing two columns both named `A1` (and `A2`). Fixed by adding `--ignore effect_allele,other_allele` (supported by this LDSC fork) to the job definition, which excludes the duplicate pair from auto-matching. The raw file had to be re-downloaded (570 MB, ~1 min) because the failed run's cleanup step deletes the raw file unconditionally. The rerun touched only this one accession.
- Munge (after fix): 1,190,321 SNPs, mean χ² = 2.068 (35.5 s)
- `--h2-cts` run time: 5 m 25.2 s (total incl. munge: ≈ 6 m 1 s; excludes the failed attempt and re-download)
- 11 tissues by `Coefficient_P_value` (ascending):

  | Rank | Tissue | P-value |
  |---|---|---|
  | 1 | BRAIN | 0.00024 |
  | 2 | ISLET | 0.63 |
  | 3 | IMMUNE | 0.67 |
  | 4 | MUSCLE | 0.74 |
  | 5 | ARTERY | 0.79 |
  | 6 | INTESTINE | 0.86 |
  | 7 | SKIN | 0.96 |
  | 8 | LUNG | 0.96 |
  | 9 | HEART | 0.99 |
  | 10 | ADIPOSE | 1.00 |
  | 11 | LIVER | 1.00 |

- **Expected tissue (BRAIN) rank: 1 of 11.**

## Summary

All four traits land their expected tissue in the top 2 of 11: three land it at
rank 1 (ulcerative colitis → IMMUNE, BMI → BRAIN, intelligence → BRAIN) and one
at rank 2 (type 1 diabetes → IMMUNE, just behind LUNG). All four munged files
clear the >200,000-SNP and mean-χ²->1 bar by a wide margin (1,190,321 SNPs
each — the full HapMap3 reference list — with mean χ² between 1.47 and 6.47).

Every accession's raw download was byte-verified against the server's
`Content-Length` and passed `gzip -t` before munging. Raw files are deleted
after a successful munge (same behavior as the original pipeline script), so
they are not retained anywhere.
