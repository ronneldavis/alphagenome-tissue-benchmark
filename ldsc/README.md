# Matched LDSC-SEG comparison

14 traits, same GWAS datasets and same 11 tissue definitions as the AlphaGenome benchmark.

- `results/*.cell_type_results.txt` — LDSC-SEG output per trait (coefficient, SE, one-tailed P)
- `scripts/` — gene-set construction from GTEx v8, annotation building, LD-score splitting,
  the download/munge pipeline, and the comparison script
- `../data/method_comparison.csv` — the joined head-to-head table

Reproduction notes: LDSC was run from the abdenlab/ldsc-python3 fork under Python 3.12 with
nine patches applied (chromosome-split reading, one-tailed P, NA handling, tab separator,
binary-mode writes, df.ix removal, gzip text mode, error handler, entry-point aliases).
Tissue LD scores were computed from 1000G Phase 3 EUR; annotations were derived from GTEx v8
median TPM with GENCODE v26lift37 coordinates, top 10% specifically expressed genes, 100kb windows.
Caveat: GTEx has no pancreatic islet, so the ISLET annotation is bulk pancreas.
