# Raw AlphaGenome score matrices

Copied on 10 Sept 2026 from the session scratchpad where `scripts/22_run_all.py` and
`scripts/11_background.py` were run, so that the analyses can be reproduced after that
temporary directory is cleaned up.

- `bg_scores.npz` — background variants (318 x 1116 histone tracks under key `hist`; RNA under `rna`)
- `T_<study>_scores.npz` / `T_<study>_variants.csv` / `T_<study>_meta.csv` — per trait-GWAS pair
- `t2d_scores.npz`, `t2d_variants.csv` — the Xue 2018 T2D replication set (scored by the earlier `08_score.py`)
- `ldl_*` — the pilot LDL run (not used in the paper's tables)

Reproduce the review-stage numbers and tables with:

    uv run --with numpy --with pandas --with scipy --with statsmodels --with matplotlib \
        python scripts/37_review_analyses.py scores .
    tectonic main.tex
