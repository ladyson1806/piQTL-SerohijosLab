# QTL Mapping on Genome

Methods-only workflow for building a genome-wide QTL overlap summary table and visualization panels.

## Inputs

Expected files in `data/`:

- `piQTL_results.csv`
- `piQTL_SNP_annotation.csv`
- `pQTL_results.csv`
- `eQTL_results.csv`
- `eQTL_SNP_annotation.csv`

## Run

From this directory:

```bash
python 01_create_qtl_summary_table.py
python 02_visualize_qtl_overlap.py --mode exact
python 02_visualize_qtl_overlap.py --mode colocal
```

## Workflow

- `01_create_qtl_summary_table.py`
  - Builds `out/tables/QTL_overlap_summary.csv` from SNP-level annotations and QTL tables.
  - Computes exact and colocal count columns for piQTL, pQTL, and eQTL.
- `02_visualize_qtl_overlap.py`
  - Generates genome-wide overlap visualizations for `exact` and `colocal` modes.
  - Supports optional output format and figure options.

## Outputs

Main output locations:

- `out/tables/QTL_overlap_summary.csv`
- `out/figures/QTL_overlap_exact.png`
- `out/figures/QTL_overlap_exact.svg`
- `out/figures/QTL_overlap_colocal.png`
- `out/figures/QTL_overlap_colocal.svg`
