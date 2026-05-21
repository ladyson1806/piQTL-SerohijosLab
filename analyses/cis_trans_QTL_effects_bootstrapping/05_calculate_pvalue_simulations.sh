#!/bin/bash

# Calculate empirical p-values for reference values against bootstrap simulations.
#
# Outputs two one-tailed tables with +1 correction:
# 1) p(sim >= ref)
# 2) p(ref < sim) = p(sim > ref)

SCRIPT_DIR=$(dirname "${BASH_SOURCE[0]}")
cd "${SCRIPT_DIR}" || exit 1

INPUT_PQTL_SIM="out/simulation_results_pQTL/simulation_summary.tsv"
INPUT_EQTL_SIM="out/simulation_results_eQTL/simulation_summary.tsv"
INPUT_SUMMARY="out/tables/cis_trans_qtl_summary_combined.csv"
OUTPUT_PVALUES_GE="out/tables/cis_trans_qtl_pvalues.csv"
OUTPUT_PVALUES_SMALLER="out/tables/cis_trans_qtl_pvalues_ref_smaller.csv"

PVALUE_SCRIPT="src/calculate_pvalue_simulations.py"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --pqtl-sim)
            INPUT_PQTL_SIM="$2"
            shift 2
            ;;
        --eqtl-sim)
            INPUT_EQTL_SIM="$2"
            shift 2
            ;;
        --summary-table)
            INPUT_SUMMARY="$2"
            shift 2
            ;;
        --output-pvalues-ge)
            OUTPUT_PVALUES_GE="$2"
            shift 2
            ;;
        --output-pvalues-smaller)
            OUTPUT_PVALUES_SMALLER="$2"
            shift 2
            ;;
        --help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --pqtl-sim PATH         pQTL simulation_summary.tsv (default: ${INPUT_PQTL_SIM})"
            echo "  --eqtl-sim PATH         eQTL simulation_summary.tsv (default: ${INPUT_EQTL_SIM})"
            echo "  --summary-table PATH    cis_trans_qtl_summary_combined.csv (default: ${INPUT_SUMMARY})"
            echo "  --output-pvalues-ge PATH       Output GE-table CSV (default: ${OUTPUT_PVALUES_GE})"
            echo "  --output-pvalues-smaller PATH  Output smaller-than-table CSV (default: ${OUTPUT_PVALUES_SMALLER})"
            echo "  --help                  Show this message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

if [ ! -f "${INPUT_PQTL_SIM}" ]; then
    echo "Error: Input file not found: ${INPUT_PQTL_SIM}"
    exit 1
fi

if [ ! -f "${INPUT_EQTL_SIM}" ]; then
    echo "Error: Input file not found: ${INPUT_EQTL_SIM}"
    exit 1
fi

if [ ! -f "${INPUT_SUMMARY}" ]; then
    echo "Error: Input file not found: ${INPUT_SUMMARY}"
    exit 1
fi

if [ ! -f "${PVALUE_SCRIPT}" ]; then
    echo "Error: Python script not found: ${PVALUE_SCRIPT}"
    exit 1
fi

echo "========================================================================"
echo "Calculate Empirical P-values for Bootstrap Simulations (+1 correction)"
echo "========================================================================"
echo ""
echo "Inputs:"
echo "  pQTL simulation: ${INPUT_PQTL_SIM}"
echo "  eQTL simulation: ${INPUT_EQTL_SIM}"
echo "  Summary table: ${INPUT_SUMMARY}"
echo ""
echo "Output:"
echo "  GE-table CSV: ${OUTPUT_PVALUES_GE}"
echo "  smaller-than-table CSV: ${OUTPUT_PVALUES_SMALLER}"
echo ""

CMD=(
    python3 "${PVALUE_SCRIPT}"
    --pqtl-sim "${INPUT_PQTL_SIM}"
    --eqtl-sim "${INPUT_EQTL_SIM}"
    --summary-table "${INPUT_SUMMARY}"
    --output-pvalues-ge "${OUTPUT_PVALUES_GE}"
    --output-pvalues-smaller "${OUTPUT_PVALUES_SMALLER}"
)

echo "Command:"
printf '  %q ' "${CMD[@]}"
echo ""
echo ""

"${CMD[@]}"
if [ $? -ne 0 ]; then
    echo "Error: P-value calculation script failed"
    exit 1
fi

echo ""
echo "Completed successfully."
echo "========================================================================"
