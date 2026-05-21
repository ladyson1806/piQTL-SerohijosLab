# piQTL-SerohijosLab

Code and reproducibility materials for QTL SNP comparison analyses performed for manuscript-associated work.

This repository is a methods and usage resource. It does not include all analyses from the manuscript, and it does not provide biological interpretation of results.

## Environment

- OS: Linux/macOS recommended.
- Runtime: Python 3.x and shell.
- Common Python packages used across directories:
	- pandas
	- numpy
	- scipy
	- matplotlib
	- seaborn
	- plotly

## QTL SNP Comparison Analyses

The following directories contain independent analysis modules for QTL SNP comparison tasks.

### 1) `/analyses/cis_trans_QTL_effects/`

Role: cis/trans effect-size processing and comparison workflows for multiple source datasets.

Contains independent dated modules:

- `2025-04-17_cis_trans_QTLs/`
- `2026-02-03_cis_trans_QTLs_in_natural_isolates/`
- `2026-02-09_cis_trans_QTLs_from_human_studies/`

### 2) `/analyses/comparing_piQTL_pQTL_eQTL/`

Role: build exact-match and LD-colocalized overlap tables across piQTL, pQTL, and eQTL studies.

Main scripts:

- `01_Jakobson_QTL_position_comparison.sh`
- `02_make_summary_piQTL_eQTL_pQTL_comparison.sh`
- `03_make_summary_piQTL_eQTL_pQTL_comparison_colocal.sh`
- `04_Add_SNP_annotation_to_summary_tables.sh`
- `05_compare_piQTL_based_on_QTL_overlaps.sh`
- `06_add_QTL_ORF_positions.sh`

### 3) `/analyses/QTL_mapping_on_genome/`

Role: generate a genome-wide summary table and visualization panels for exact/colocal QTL overlap counts.

Main scripts:

- `01_main.sh`
- `01_create_qtl_summary_table.py`
- `02_visualize_qtl_overlap.py`

### 4) `/analyses/QTL_SNPs_overlap_enrichments_bootstrapping/`

Role: test overlap enrichment with randomized SNP sampling and empirical p-value calculation.

Main scripts:

- `01_standardize_tables.sh`
- `02_randomize_overlap_counts.sh`
- `03_make_SNP_status_table.sh`
- `04_make_SNP_overlap_count_histgrams.sh`
- `05_calculate_empirical_pvalues.sh`

Alternative scripts for shared-SNP-specific analysis are also provided in this directory (`02b`, `03b`, `04b`, `05b`).

### 5) `cis_trans_QTL_effects_bootstrapping/`

Role: bootstrap simulation and summary-statistics comparison of cis/trans effect-size metrics.

Main scripts:

- `01a_random_sampling_pQTL_simulation.sh`
- `01b_random_sampling_eQTL_simulation.sh`
- `02_filter_and_plot_piQTL_based_on_cis_condition.sh`
- `03_summarize_cis_trans_qtl_effects.sh`
- `04_compare_pqtl_eqtl_simulations.sh`
- `05_calculate_pvalue_simulations.sh`

### Notes

- Each directory contains a local README with concise usage for that module.
- Input data are expected in each module's `data/` directory.
