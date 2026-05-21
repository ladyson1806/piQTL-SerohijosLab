# cis/trans QTL Effects Bootstrapping

Methods-only workflow for bootstrap simulation and comparative summary metrics across piQTL, pQTL, and eQTL cis/trans effect-size datasets.

## Inputs

Expected files in `data/`:

- `Jakobson_tested_genes.csv`
- `Jakobson_pQTLs.tsv`
- `Albert_eQTLs.tsv`
- `piQTL_results.csv`
- `piQTL_target_genes.csv`

## Run

From this directory:

```bash
bash 01a_random_sampling_pQTL_simulation.sh
bash 01b_random_sampling_eQTL_simulation.sh
bash 02_filter_and_plot_piQTL_based_on_cis_condition.sh
bash 03_summarize_cis_trans_qtl_effects.sh
bash 04_compare_pqtl_eqtl_simulations.sh
bash 05_calculate_pvalue_simulations.sh
```

## Workflow

- `01a_random_sampling_pQTL_simulation.sh`
  - Runs pQTL bootstrap random sampling.
  - Generates simulation summary tables and figure outputs.
- `01b_random_sampling_eQTL_simulation.sh`
  - Runs eQTL bootstrap random sampling with analogous outputs.
- `02_filter_and_plot_piQTL_based_on_cis_condition.sh`
  - Generates piQTL-side filtered summaries and figure assets used in comparison steps.
- `03_summarize_cis_trans_qtl_effects.sh`
  - Produces global and piQTL-target subset summary metrics for piQTL, pQTL, and eQTL.
- `04_compare_pqtl_eqtl_simulations.sh`
  - Builds panel plots comparing simulation distributions against observed references.
- `05_calculate_pvalue_simulations.sh`
  - Calculates empirical p-values for simulation-versus-reference comparisons.

## Key Metric Definitions

Computed metrics include:

- `mean_diff = mean(cis) - mean(trans)`
- `median_diff = median(cis) - median(trans)`
- `cis_trans_ratio = n_cis / n_trans`
- `trans_pct = 100 * n_trans / (n_cis + n_trans)`

Empirical p-value formula:

- `p = (n_extreme + 1) / (n_total + 1)`

## Outputs

Main output locations:

- `out/simulation_results_pQTL/`
- `out/simulation_results_eQTL/`
- `out/tables/`
- `out/simulation_comparison/`

Representative outputs:

- `out/tables/cis_trans_qtl_summary_global.csv`
- `out/tables/cis_trans_qtl_summary_piqtl_target_subset.csv`
- `out/tables/cis_trans_qtl_summary_combined.csv`
- `out/tables/cis_trans_qtl_pvalues.csv`
- `out/tables/cis_trans_qtl_pvalues_ref_smaller.csv`
