#!/bin/bash

# Compare pQTL and eQTL bootstrap simulation results.
#
# Creates 4 figure panels showing distributions from 10,000 random samplings
# with reference markers for global, piQTL-target subset, and piQTL values.
# Output: PNG and SVG formats only

SCRIPT_DIR=$(dirname "${BASH_SOURCE[0]}")
cd "${SCRIPT_DIR}" || exit 1

INPUT_PQTL_SIM="out/simulation_results_pQTL/simulation_summary.tsv"
INPUT_EQTL_SIM="out/simulation_results_eQTL/simulation_summary.tsv"
INPUT_SUMMARY="out/tables/cis_trans_qtl_summary_combined.csv"

COMPARISON_SCRIPT="src/compare_pqtl_eqtl_simulations.py"
OUTPUT_DIR="out/simulation_comparison"

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
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --pqtl-sim PATH         pQTL simulation_summary.tsv (default: ${INPUT_PQTL_SIM})"
            echo "  --eqtl-sim PATH         eQTL simulation_summary.tsv (default: ${INPUT_EQTL_SIM})"
            echo "  --summary-table PATH    cis_trans_qtl_summary_combined.csv (default: ${INPUT_SUMMARY})"
            echo "  --output-dir PATH       Output directory (default: ${OUTPUT_DIR})"
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

if [ ! -f "${COMPARISON_SCRIPT}" ]; then
    echo "Error: Python script not found: ${COMPARISON_SCRIPT}"
    exit 1
fi

echo "========================================================================"
echo "Compare pQTL and eQTL Simulation Results"
echo "========================================================================"
echo ""
echo "Inputs:"
echo "  pQTL simulation: ${INPUT_PQTL_SIM}"
echo "  eQTL simulation: ${INPUT_EQTL_SIM}"
echo "  Summary table: ${INPUT_SUMMARY}"
echo ""
echo "Output:"
echo "  Directory: ${OUTPUT_DIR}"
echo "  Formats: PNG, SVG"
echo ""

CMD=(
    python3 "${COMPARISON_SCRIPT}"
    --pqtl-sim "${INPUT_PQTL_SIM}"
    --eqtl-sim "${INPUT_EQTL_SIM}"
    --summary-table "${INPUT_SUMMARY}"
    --output-dir "${OUTPUT_DIR}"
)

echo "Command:"
printf '  %q ' "${CMD[@]}"
echo ""
echo ""

"${CMD[@]}"
if [ $? -ne 0 ]; then
    echo "Error: Comparison script failed"
    exit 1
fi

echo ""
echo "Completed successfully."
echo "========================================================================"
