# cis/trans QTL Effects (2025-04-17)

Baseline workflow to identify cis/trans classes and generate effect-size comparison outputs for multiple QTL datasets.

## Inputs

Expected files in `data/`:

- `2025-04-17_piQTL_SNPs.csv`
- `2025-04-17_snps_annotations_LD050_cisSNP-5kb.csv`
- `2025-04-17_PPI_gene_list.csv`
- `2025-04-17_Chick_pQTL_Table_S5.csv`
- `2025-04-18_Chick_eQTL_Table_S4.csv`
- `Jakobson_pQTLs.csv` (for Jakobson-specific plotting step)

## Run

From this directory:

```bash
bash 01_main.sh
bash 02_Jakobson_boxplot.sh
```

## What the scripts do

- `01_main.sh`
  - Generates formatted cis/trans tables for piQTL, Albert eQTL, Chick pQTL, and Chick eQTL.
  - Produces boxplots and t-test result files for these datasets.
- `02_Jakobson_boxplot.sh`
  - Generates formatted cis/trans Jakobson pQTL tables.
  - Produces boxplots and t-test outputs for all genes and piQTL target-gene subsets.

## Outputs

Main output locations:

- `out/formatted_tables/`
- `out/figures/`
- `out/t_test/`
