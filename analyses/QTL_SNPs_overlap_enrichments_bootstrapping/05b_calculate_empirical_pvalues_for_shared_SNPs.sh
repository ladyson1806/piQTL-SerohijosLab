#!/bin/sh

cd $(dirname ${0})

# Calculate empirical p-values for QTL overlap enrichment

mkdir -p out/empirical_pvalues_for_shared_snps

RANDOMIZED_EXACT="out/randomized_overlap_counts_exact_shared_snps/overlap_counts_results.csv"
RANDOMIZED_COLOCAL="out/randomized_overlap_counts_colocal_shared_snps/overlap_counts_results.csv"
ACTUAL_CSV="data/actual_piQTL_SNPs_overlap_status_for_shared_SNPs.csv"
OUTPUT="out/empirical_pvalues_for_shared_snps/empirical_pvalues_summary.csv"

python src/calculate_empirical_pvalues.py \
	--randomized-exact ${RANDOMIZED_EXACT} \
	--randomized-colocal ${RANDOMIZED_COLOCAL} \
	--actual ${ACTUAL_CSV} \
	--output ${OUTPUT}

echo ""
echo "Done! Empirical p-values saved to: ${OUTPUT}"
