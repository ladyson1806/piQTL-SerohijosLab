# cis/trans QTL Effects from Natural Isolates

Methods-only workflow for generating cis/trans-formatted tables, boxplots, and t-test outputs from natural-isolate QTL datasets.

## Inputs

Expected files in `data/`:

- `Teyssonniere_pQTLs.csv`
- `Teyssonniere_eQTLs.csv`
- `Loegler_QTLs.csv`
- `PPI_gene_list.csv`

## Run

From this directory:

```bash
bash 01_boxplot_for_Teyssonniere.sh
bash 02_boxplot_for_Loegler.sh
```

## Workflow

- `01_boxplot_for_Teyssonniere.sh`
  - Formats Teyssonniere pQTL and eQTL tables with cis/trans labels.
  - Produces all-genes and piQTL-target-gene subset outputs.
  - Generates boxplots and t-test result files.
- `02_boxplot_for_Loegler.sh`
  - Extracts and formats Proteomics and Transcriptomics rows from Loegler data.
  - Produces all-genes and piQTL-target-gene subset outputs.
  - Generates boxplots and t-test result files.

## Outputs

Main output locations:

- `out/formatted_tables/`
- `out/figures/`
- `out/t_test/`

Representative files:

- `out/formatted_tables/Teyssonniere_pQTLs_formatted.csv`
- `out/formatted_tables/Loegler_eQTLs_formatted.csv`
- `out/figures/boxplot_Teyssonniere_pQTLs.png`
- `out/figures/boxplot_Loegler_eQTLs.png`
