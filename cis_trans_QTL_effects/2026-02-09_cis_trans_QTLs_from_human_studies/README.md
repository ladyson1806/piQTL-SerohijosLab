# cis/trans QTL Effects from Human Studies

Methods-only workflow for generating cis/trans-formatted tables, boxplots, and t-test outputs from human QTL studies.

## Inputs

Expected files in `data/`:

- `Li_pQTLs.csv`
- `Vosa_cis_eQTLs.tsv.gz`
- `Vosa_trans_eQTLs.tsv`

## Run

From this directory:

```bash
bash 01_preprocess.sh
bash 02_boxplot_for_Li.sh
bash 03_boxplot_for_Vosa.sh
bash 04_boxplot_for_Li_by_Zscore.sh
```

## Workflow

- `01_preprocess.sh`
  - Creates required output directories.
  - Expands `Vosa_cis_eQTLs.tsv.gz` into `data/Vosa_cis_eQTLs.tsv`.
- `02_boxplot_for_Li.sh`
  - Formats Li pQTL data.
  - Generates absolute-effect boxplot and t-test output.
- `03_boxplot_for_Vosa.sh`
  - Merges Vosa cis/trans files into a single formatted eQTL table.
  - Generates standard and zoomed boxplots plus t-test output.
- `04_boxplot_for_Li_by_Zscore.sh`
  - Runs Li analysis using Z-score-based formatting.
  - Generates Z-score boxplot and t-test output.

## Outputs

Main output locations:

- `out/formatted_tables/`
- `out/figures/`
- `out/t_test/`

Representative files:

- `out/formatted_tables/Li_pQTLs_formatted.csv`
- `out/formatted_tables/Vosa_eQTLs_formatted.csv`
- `out/figures/boxplot_Li_pQTLs.png`
- `out/figures/boxplot_Vosa_eQTLs.png`
