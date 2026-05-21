#!/bin/bash

# Random sampling simulation for cis and trans eQTL effect distribution analysis
#
# This script runs Monte Carlo simulations to analyze how cis and trans eQTL
# effect sizes differ when randomly sampling genes from the full gene universe.
#
# Output:
# - simulation_summary.tsv: 10,000 rows × 9 columns (metrics for each simulation)
# - gene_samples.tsv: 10,000 rows × 44 columns (genes sampled in each simulation)
# - boxplot_figure.html/png/svg: 5-panel visualization
# - checkpoints/boxplot_figure_cp*.png/svg: grouped cis/trans effect size boxplots every 100 simulations
# - summary_statistics.txt: Descriptive statistics for all metrics

# Get the directory where this script is located
SCRIPT_DIR=$(dirname "${BASH_SOURCE[0]}")
cd "${SCRIPT_DIR}" || exit 1

# Define file paths
INPUT_EQTL="data/Albert_eQTLs.tsv"
OUTPUT_DIR="out/simulation_results_eQTL"
SIMULATION_SCRIPT="src/random_sampling_eQTL_simulation.py"

# Check if input files exist
if [ ! -f "${INPUT_EQTL}" ]; then
    echo "Error: Input file not found: ${INPUT_EQTL}"
    exit 1
fi

if [ ! -f "${SIMULATION_SCRIPT}" ]; then
    echo "Error: Simulation script not found: ${SIMULATION_SCRIPT}"
    exit 1
fi

# Parse command-line arguments
N_GENES=44
N_SIMULATIONS=10000
RANDOM_SEED=42
CUSTOM_OUTPUT=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --n-genes)
            N_GENES=$2
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
            echo "  --n-genes NUM                 Number of genes to sample per simulation (default: 44)"
            echo "  --n-simulations NUM           Number of simulations to run (default: 10000)"
            echo "  --random-seed NUM             Random seed for reproducibility (default: 42)"
            echo "  --output-dir PATH             Output directory (default: out/simulation_results_eQTL)"
            echo "  --help                        Show this help message"
            echo ""
            echo "Examples:"
            echo "  # Default run"
            echo "  $0"
            echo ""
            echo "  # Quick test with small sample"
            echo "  $0 --n-genes 10 --n-simulations 10 --output-dir out/simulation_results_eQTL/test"
            echo ""
            echo "  # Custom parameters with different seed"
            echo "  $0 --n-genes 50 --n-simulations 500 --random-seed 123"
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
echo "Random Sampling eQTL Simulation"
echo "========================================================================"
echo ""
echo "Parameters:"
echo "  Number of genes per simulation: ${N_GENES}"
echo "  Number of simulations: ${N_SIMULATIONS}"
echo "  Random seed: ${RANDOM_SEED}"
echo "  Output directory: ${OUTPUT_DIR}"
echo ""
echo "Input files:"
echo "  eQTL data: ${INPUT_EQTL}"
echo ""

# Run the simulation
echo "Starting simulation..."
python "${SIMULATION_SCRIPT}" \
    --eqtl-file "${INPUT_EQTL}" \
    --n-genes "${N_GENES}" \
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
echo "  - gene_samples.tsv"
echo "    ${N_SIMULATIONS} simulations × ${N_GENES} genes per simulation"
echo "    Track which genes were sampled in each simulation"
echo ""
echo "  - boxplot_figure (html/png/svg)"
echo "    5-panel boxplot visualization:"
echo "      1. Mean effect size difference (cis - trans)"
echo "      2. Median effect size difference (cis - trans)"
echo "      3. cis-eQTL count distribution"
echo "      4. trans-eQTL count distribution"
echo "      5. cis/trans eQTL ratio distribution"
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
