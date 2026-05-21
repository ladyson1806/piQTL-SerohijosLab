# QTL SNP Overlap Enrichment Bootstrapping

Methods-only workflow for standardizing QTL tables, randomizing overlap counts, and calculating empirical p-values.

## Overlap Definitions

- Exact overlap: piQTL `SNP_marker` matches pQTL/eQTL `SNP_marker`.
- Colocal overlap: piQTL and QTL peaks are mutually contained in each other's interval ranges on the same chromosome.

## Inputs

- `data/actual_piQTL_SNPs.txt`
- `data/masked_SNPs.txt` (used during standardization)
- raw tables required by the standardization scripts

After `01_standardize_tables.sh`, downstream scripts use:

- `out/standardized_tables/piQTL_SNP_annotation.csv`
- `out/standardized_tables/pQTL_results.csv`
- `out/standardized_tables/eQTL_results.csv`

## Run

From this directory:

```bash
bash 01_standardize_tables.sh
bash 02_randomize_overlap_counts.sh
bash 03_make_SNP_status_table.sh
bash 04_make_SNP_overlap_count_histgrams.sh
bash 05_calculate_empirical_pvalues.sh
```

Shared-SNP variants are provided as alternatives:

- `02b_randomize_overlap_counts_only_for_shared_SNPs.sh`
- `03b_make_SNP_status_table_for_shared_SNPs.sh`
- `04b_make_SNP_overlap_count_histgrams_shared_SNPs.sh`
- `05b_calculate_empirical_pvalues_for_shared_SNPs.sh`

## Workflow

- `01_standardize_tables.sh`
  - Standardizes piQTL, pQTL, and eQTL tables into a common schema.
- `02_randomize_overlap_counts.sh`
  - Runs random sampling for exact/colocal overlap counts.
- `03_make_SNP_status_table.sh`
  - Builds per-SNP overlap status tables.
- `04_make_SNP_overlap_count_histgrams.sh`
  - Generates histogram panels for randomized overlap distributions.
- `05_calculate_empirical_pvalues.sh`
  - Calculates enrichment p-values and fold-enrichment summaries.

Empirical p-value formula:

`p = (n_extreme + 1) / (n_total + 1)`

## Outputs

Main output locations:

- `out/standardized_tables/`
- `out/randomized_overlap_counts_exact/`
- `out/randomized_overlap_counts_colocal/`
- `out/overlap_status/`
- `out/overlap_histograms/`
- `out/empirical_pvalues/`
