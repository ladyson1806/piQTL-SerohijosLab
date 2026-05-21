# cis_trans_QTL_effects

Independent cis/trans QTL effect-size workflows grouped by analysis date and source dataset.

## Directory Map

- `2025-04-17_cis_trans_QTLs/`
  - Baseline cis/trans processing for piQTL, Albert eQTL, Chick pQTL/eQTL, and Jakobson pQTL follow-up plotting.
- `2026-02-03_cis_trans_QTLs_in_natural_isolates/`
  - cis/trans effect-size workflow for Teyssonniere and Loegler natural-isolate datasets.
- `2026-02-09_cis_trans_QTLs_from_human_studies/`
  - cis/trans effect-size workflow for human studies (Li pQTL and Vosa eQTL).

## Usage

Use each dated subdirectory as an independent module:

1. Move into the target subdirectory.
2. Ensure required input files are available in `data/`.
3. Run the local shell scripts listed in that subdirectory README.
4. Check outputs in `out/formatted_tables/`, `out/figures/`, and `out/t_test/`.

## Notes

- This directory does not define a single global run order across dated modules.
- cis/trans definitions are implemented in each module's formatting scripts and reflect source-specific table conventions.
