#!/bin/bash

# Summarize cis/trans QTL effect-size metrics for piQTL, pQTL, and eQTL.
#
# Produces:
# - Global summary (all rows)
# - piQTL-target subset summary
# - Optional combined summary table

SCRIPT_DIR=$(dirname "${BASH_SOURCE[0]}")
cd "${SCRIPT_DIR}" || exit 1

INPUT_PIQTL="out/tables/piQTL_results_cis_trans_annotated.csv"
INPUT_PQTL="data/Jakobson_pQTLs.tsv"
INPUT_EQTL="data/Albert_eQTLs.tsv"
INPUT_TARGET="data/piQTL_target_genes.csv"

SUMMARY_SCRIPT="src/summarize_cis_trans_qtl_effects.py"

OUT_GLOBAL="out/tables/cis_trans_qtl_summary_global.csv"
OUT_TARGET="out/tables/cis_trans_qtl_summary_piqtl_target_subset.csv"
OUT_COMBINED="out/tables/cis_trans_qtl_summary_combined.csv"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --piqtl-file)
            INPUT_PIQTL="$2"
            shift 2
            ;;
        --pqtl-file)
            INPUT_PQTL="$2"
            shift 2
            ;;
        --eqtl-file)
            INPUT_EQTL="$2"
            shift 2
            ;;
        --target-file)
            INPUT_TARGET="$2"
            shift 2
            ;;
        --out-global)
            OUT_GLOBAL="$2"
            shift 2
            ;;
        --out-target)
            OUT_TARGET="$2"
            shift 2
            ;;
        --out-combined)
            OUT_COMBINED="$2"
            shift 2
            ;;
        --help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --piqtl-file PATH     Annotated piQTL CSV (default: ${INPUT_PIQTL})"
            echo "  --pqtl-file PATH      pQTL TSV (default: ${INPUT_PQTL})"
            echo "  --eqtl-file PATH      eQTL TSV (default: ${INPUT_EQTL})"
            echo "  --target-file PATH    piQTL target gene CSV (default: ${INPUT_TARGET})"
            echo "  --out-global PATH     Output CSV for global summary (default: ${OUT_GLOBAL})"
            echo "  --out-target PATH     Output CSV for target subset summary (default: ${OUT_TARGET})"
            echo "  --out-combined PATH   Output CSV for combined summary (default: ${OUT_COMBINED})"
            echo "  --help                Show this message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

if [ ! -f "${INPUT_PIQTL}" ]; then
    echo "Error: Input file not found: ${INPUT_PIQTL}"
    exit 1
fi

if [ ! -f "${INPUT_PQTL}" ]; then
    echo "Error: Input file not found: ${INPUT_PQTL}"
    exit 1
fi

if [ ! -f "${INPUT_EQTL}" ]; then
    echo "Error: Input file not found: ${INPUT_EQTL}"
    exit 1
fi

if [ ! -f "${INPUT_TARGET}" ]; then
    echo "Error: Input file not found: ${INPUT_TARGET}"
    exit 1
fi

if [ ! -f "${SUMMARY_SCRIPT}" ]; then
    echo "Error: Python script not found: ${SUMMARY_SCRIPT}"
    exit 1
fi

echo "========================================================================"
echo "Cis/Trans QTL Effect Summary"
echo "========================================================================"
echo ""
echo "Inputs:"
echo "  piQTL (annotated): ${INPUT_PIQTL}"
echo "  pQTL: ${INPUT_PQTL}"
echo "  eQTL: ${INPUT_EQTL}"
echo "  piQTL target genes: ${INPUT_TARGET}"
echo ""
echo "Outputs:"
echo "  Global summary: ${OUT_GLOBAL}"
echo "  piQTL-target subset summary: ${OUT_TARGET}"
echo "  Combined summary: ${OUT_COMBINED}"
echo ""

CMD=(
    python3 "${SUMMARY_SCRIPT}"
    --piqtl-file "${INPUT_PIQTL}"
    --pqtl-file "${INPUT_PQTL}"
    --eqtl-file "${INPUT_EQTL}"
    --target-file "${INPUT_TARGET}"
    --out-global "${OUT_GLOBAL}"
    --out-target "${OUT_TARGET}"
    --out-combined "${OUT_COMBINED}"
)

echo "Command:"
printf '  %q ' "${CMD[@]}"
echo ""
echo ""

"${CMD[@]}"
if [ $? -ne 0 ]; then
    echo "Error: Summary script failed"
    exit 1
fi

echo ""
echo "Completed successfully."
echo "========================================================================"
