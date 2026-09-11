# Per-variant tissue assignment by a sequence-to-function model is limited by effect detectability more than by assay coverage

Author: Ronnel Davis

Benchmarks AlphaGenome per-variant tissue assignments against 20 trait-GWAS
pairs with pre-specified expected tissues. Includes a matched LDSC-SEG
comparison run on the same traits and tissue definitions.

## Directory map

- `main.tex` — manuscript source (compiles to `main.pdf`)
- `data/` — benchmark results and metadata, CSV/JSON (~0.5 MB)
- `scores/` — raw AlphaGenome score matrices (~13 MB); see `scores/README.md`
- `scripts/` — Python analysis scripts
- `ldsc/` — LDSC-SEG comparison results, scripts, and README
- `REVIEW.md` — internal review notes

## Reproducing the results

From the repository root:

    uv run --with numpy --with pandas --with scipy --with statsmodels --with matplotlib \
        python scripts/37_review_analyses.py scores .
    tectonic main.tex

See `scores/README.md` for what each file in `scores/` contains.

## LDSC-SEG comparison

The LDSC-SEG reference data and the patched LDSC fork used to run it are not
included in this repository. See `ldsc/README.md` for what was used and how
it was run.

## Data availability

All GWAS inputs are public; see Appendix Table A of the manuscript for accessions.

## License

License: to be decided by the author.
