#!/bin/sh

cd $(dirname ${0})

mkdir -p out/piQTL_effect_sizes
mkdir -p out/piQTL_effect_sizes/figures


# Define paths to input files
PIQTL_EFFECT_SIZES="data/piQTLs_formatted_lead.csv"
PIQTL_OVERLAPS_EXACT="out/summary/summary_master_SNP_tabl.csv"
PIQTL_OVERLAPS_COLOCAL="out/summary_colocal/summary_master_SNP_table_grouped_by_piQTL_full.csv"
OUTPUT_DIR="out/piQTL_effect_sizes/figures"

# Compare piQTL effect sizes based on exact QTL overlaps
EXACT_OVERLAP_SUMMARY_TABLE="out/piQTL_effect_sizes/exact_overlap_summary_table.csv"
EXACT_OVERLAP_FIGURE="out/piQTL_effect_sizes/figures/exact_overlap_effect_size_comparison.png"
EXACT_OVERLAP_STATISTICS="out/piQTL_effect_sizes/exact_overlap_effect_size_statistics.txt"
if [ ! -f "$EXACT_OVERLAP_SUMMARY_TABLE" ] || [ ! -f "$EXACT_OVERLAP_FIGURE" ]; then
  python src/classify_and_plot_piQTL_overlaps.py \
    --piqtl_input "$PIQTL_EFFECT_SIZES" \
    --overlap_input "$PIQTL_OVERLAPS_EXACT" \
    --output_table "$EXACT_OVERLAP_SUMMARY_TABLE" \
    --output_figure "$EXACT_OVERLAP_FIGURE" \
    --output_stats "$EXACT_OVERLAP_STATISTICS"
fi


# Compare piQTL effect sizes based on colocalized QTL overlaps
COLOCAL_OVERLAP_SUMMARY_TABLE="out/piQTL_effect_sizes/colocal_overlap_summary_table.csv"
COLOCAL_OVERLAP_FIGURE="out/piQTL_effect_sizes/figures/colocal_overlap_effect_size_comparison.png"
COLOCAL_OVERLAP_STATISTICS="out/piQTL_effect_sizes/colocal_overlap_effect_size_statistics.txt"
if [ ! -f "$COLOCAL_OVERLAP_SUMMARY_TABLE" ] || [ ! -f "$COLOCAL_OVERLAP_FIGURE" ]; then
  python src/classify_and_plot_piQTL_overlaps.py \
    --piqtl_input "$PIQTL_EFFECT_SIZES" \
    --overlap_input "$PIQTL_OVERLAPS_COLOCAL" \
    --output_table "$COLOCAL_OVERLAP_SUMMARY_TABLE" \
    --output_figure "$COLOCAL_OVERLAP_FIGURE" \
    --output_stats "$COLOCAL_OVERLAP_STATISTICS"
fi
