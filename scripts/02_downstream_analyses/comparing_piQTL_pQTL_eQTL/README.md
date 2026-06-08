# Comparing piQTL, pQTL, and eQTL

Methods-only workflow for generating exact-match and colocalized overlap tables across piQTL, pQTL, and eQTL datasets.

## Overlap Definitions

- Exact overlap: lead SNP positions are identical across datasets.
- Colocal overlap: lead SNPs are treated as overlapping if they fall within matched LD-supported interval definitions used in this workflow.

## Inputs

Expected files are organized under `data/` and include:

- piQTL formatted lead SNP table
- Jakobson pQTL data and LD metadata
- Albert eQTL formatted table
- annotation resources used by the summary-annotation step

## Run

From this directory:

```bash
bash 01_Jakobson_QTL_position_comparison.sh
bash 02_make_summary_piQTL_eQTL_pQTL_comparison.sh
bash 03_make_summary_piQTL_eQTL_pQTL_comparison_colocal.sh
bash 04_Add_SNP_annotation_to_summary_tables.sh
bash 05_compare_piQTL_based_on_QTL_overlaps.sh
bash 06_add_QTL_ORF_positions.sh
```

## Workflow

- `01_Jakobson_QTL_position_comparison.sh`
  - Parses Jakobson LD block ranges.
  - Adds LD block information to pQTL rows.
  - Generates pairwise and three-way overlap tables.
  - Builds directly affected-gene overlap tables for exact and colocal modes.
- `02_make_summary_piQTL_eQTL_pQTL_comparison.sh`
  - Builds the master SNP summary table in exact mode.
- `03_make_summary_piQTL_eQTL_pQTL_comparison_colocal.sh`
  - Builds grouped colocal summary tables.
- `04_Add_SNP_annotation_to_summary_tables.sh`
  - Adds annotation columns to summary outputs.
- `05_compare_piQTL_based_on_QTL_overlaps.sh` and `06_add_QTL_ORF_positions.sh`
  - Provide additional comparison/annotation transformations used in downstream tables.

## Outputs

Main output locations:

- `out/piQTL_vs_pQTL/`
- `out/pQTL_vs_eQTL/`
- `out/piQTLs_vs_eQTLs/`
- `out/piQTL_vs_eQTL_vs_pQTL/`
- `out/summary/`
- `out/summary_colocal/`
- `out/summary_annotated/`
- `out/formatted_tables/`
- `out/t_test/`
