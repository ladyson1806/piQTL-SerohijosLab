#!/bin/sh

cd $(dirname ${0})

mkdir -p out/overlap_histograms_shared_snps

OUTPUT_DIR="out/overlap_histograms_shared_snps"

# Generate histogram panels for QTL overlap randomized simulation results
# Each figure panel shows 3 subplots for pQTL, eQTL, and both overlaps

RANDOMIZED_EXACT="out/randomized_overlap_counts_exact_shared_snps/overlap_counts_results.csv"
RANDOMIZED_COLOCAL="out/randomized_overlap_counts_colocal_shared_snps/overlap_counts_results.csv"
ACTUAL_CSV="data/actual_piQTL_SNPs_overlap_status_for_shared_SNPs.csv"

# Generate figure panel for EXACT overlap
python src/plot_overlap_histograms.py \
	--overlap-mode exact \
	--randomized-csv ${RANDOMIZED_EXACT} \
	--actual-csv ${ACTUAL_CSV} \
	--output-prefix ${OUTPUT_DIR}/overlap_histograms_exact

# Generate figure panel for COLOCAL overlap
python src/plot_overlap_histograms.py \
	--overlap-mode colocal \
	--randomized-csv ${RANDOMIZED_COLOCAL} \
	--actual-csv ${ACTUAL_CSV} \
	--output-prefix ${OUTPUT_DIR}/overlap_histograms_colocal

echo "Done! Generated histograms:"
echo "  - ${OUTPUT_DIR}/overlap_histograms_exact.png"
echo "  - ${OUTPUT_DIR}/overlap_histograms_exact.svg"
echo "  - ${OUTPUT_DIR}/overlap_histograms_colocal.png"
echo "  - ${OUTPUT_DIR}/overlap_histograms_colocal.svg"
