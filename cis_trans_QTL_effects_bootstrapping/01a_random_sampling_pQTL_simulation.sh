#!/bin/bash

# Random sampling simulation for cis and trans pQTL effect distribution analysis
#
# This script runs Monte Carlo simulations to analyze how cis and trans pQTL
# effect sizes differ when randomly sampling proteins from the full protein universe.
#
# Output:
# - simulation_summary.tsv: 10,000 rows × 9 columns (metrics for each simulation)
# - protein_samples.tsv: 10,000 rows × 44 columns (proteins sampled in each simulation)
# - boxplot_figure.html/png/svg: 5-panel visualization
# - checkpoints/boxplot_figure_cp*.png/svg: grouped cis/trans effect size boxplots every 100 simulations
# - summary_statistics.txt: Descriptive statistics for all metrics

# Get the directory where this script is located
SCRIPT_DIR=$(dirname "${BASH_SOURCE[0]}")
cd "${SCRIPT_DIR}" || exit 1

# Define file paths
INPUT_PROTEINS="data/Jakobson_tested_genes.csv"
INPUT_PQTL="data/Jakobson_pQTLs.tsv"
OUTPUT_DIR="out/simulation_results_pQTL"
SIMULATION_SCRIPT="src/random_sampling_pQTL_simulation.py"

# Check if input files exist
if [ ! -f "${INPUT_PROTEINS}" ]; then
    echo "Error: Input file not found: ${INPUT_PROTEINS}"
    exit 1
fi

if [ ! -f "${INPUT_PQTL}" ]; then
    echo "Error: Input file not found: ${INPUT_PQTL}"
    exit 1
fi

if [ ! -f "${SIMULATION_SCRIPT}" ]; then
    echo "Error: Simulation script not found: ${SIMULATION_SCRIPT}"
    exit 1
fi

# Parse command-line arguments
N_PROTEINS=44
N_SIMULATIONS=10000
RANDOM_SEED=42
CUSTOM_OUTPUT=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --n-proteins)
            N_PROTEINS=$2
            shift 2
            ;;
        --n-simulations)
            N_SIMULATIONS=$2
            shift 2
            ;;
        --random-seed)
            RANDOM_SEED=$2
            shift 2
            ;;
        --output-dir)
            CUSTOM_OUTPUT=$2
            shift 2
            ;;
        --help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --n-proteins NUM              Number of proteins to sample per simulation (default: 44)"
            echo "  --n-simulations NUM           Number of simulations to run (default: 10000)"
            echo "  --random-seed NUM             Random seed for reproducibility (default: 42)"
            echo "  --output-dir PATH             Output directory (default: simulation_results_pQTL)"
            echo "  --help                        Show this help message"
            echo ""
            echo "Examples:"
            echo "  # Default run"
            echo "  $0"
            echo ""
            echo "  # Quick test with small sample"
            echo "  $0 --n-proteins 10 --n-simulations 10 --output-dir simulation_results/test"
            echo ""
            echo "  # Custom parameters with different seed"
            echo "  $0 --n-proteins 50 --n-simulations 500 --random-seed 123"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Use custom output directory if specified, otherwise use default
if [ -z "${CUSTOM_OUTPUT}" ]; then
    OUTPUT_DIR="${OUTPUT_DIR}"
else
    OUTPUT_DIR="${CUSTOM_OUTPUT}"
fi

# Print summary of parameters
echo "========================================================================"
echo "Random Sampling pQTL Simulation"
echo "========================================================================"
echo ""
echo "Parameters:"
echo "  Number of proteins per simulation: ${N_PROTEINS}"
echo "  Number of simulations: ${N_SIMULATIONS}"
echo "  Random seed: ${RANDOM_SEED}"
echo "  Output directory: ${OUTPUT_DIR}"
echo ""
echo "Input files:"
echo "  Protein list: ${INPUT_PROTEINS}"
echo "  pQTL data: ${INPUT_PQTL}"
echo ""

# Run the simulation
echo "Starting simulation..."
python "${SIMULATION_SCRIPT}" \
    --protein-file "${INPUT_PROTEINS}" \
    --pqtl-file "${INPUT_PQTL}" \
    --n-proteins "${N_PROTEINS}" \
    --n-simulations "${N_SIMULATIONS}" \
    --random-seed "${RANDOM_SEED}" \
    --output-dir "${OUTPUT_DIR}"

if [ $? -ne 0 ]; then
    echo "Error: Simulation failed"
    exit 1
fi

echo ""
echo "========================================================================"
echo "Simulation completed successfully!"
echo "========================================================================"
echo ""
echo "Output files generated in: ${OUTPUT_DIR}/"
echo ""
echo "Results:"
echo "  - simulation_summary.tsv"
echo "    ${N_SIMULATIONS} simulations × 9 metrics per simulation"
echo "    Includes: mean/median effect sizes, counts, differences, ratios"
echo ""
echo "  - protein_samples.tsv"
echo "    ${N_SIMULATIONS} simulations × ${N_PROTEINS} proteins per simulation"
echo "    Track which proteins were sampled in each simulation"
echo ""
echo "  - boxplot_figure (html/png/svg)"
echo "    5-panel boxplot visualization:"
echo "      1. Mean effect size difference (cis - trans)"
echo "      2. Median effect size difference (cis - trans)"
echo "      3. cis-pQTL count distribution"
echo "      4. trans-pQTL count distribution"
echo "      5. cis/trans pQTL ratio distribution"
echo ""
echo "  - summary_statistics.txt"
echo "    Descriptive statistics (mean, median, SD, min, max, quartiles)"
echo "    for each of the 5 key metrics"
echo ""
echo "  - checkpoints/boxplot_figure_cp*.png/svg"
echo "    Grouped single-panel cis/trans effect size boxplots"
echo "    Generated every 100 simulations (cp100, cp200, ..., cp1000)"
echo ""
echo "Next steps:"
echo "  1. View interactive plot: open ${OUTPUT_DIR}/boxplot_figure.html in browser"
echo "  2. Quick-check checkpoints: ls ${OUTPUT_DIR}/checkpoints/boxplot_figure_cp*.png"
echo "  3. Analyze results: load ${OUTPUT_DIR}/simulation_summary.tsv in R/Python"
echo "  4. Review statistics: cat ${OUTPUT_DIR}/summary_statistics.txt"
echo ""
echo "========================================================================"
