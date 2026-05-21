#!/bin/bash

################################################################################
# Filter and plot piQTL data based on cis-condition
#
# This script filters piQTL results to keep only PPIs with at least one
# cis-piQTL, then creates visualizations and statistical analysis.
#
# Usage: bash 03_filter_and_plot_piQTL_based_on_cis_condition.sh
#
# Author: Analysis Pipeline
# Date: 2026-02-19
################################################################################

set -e  # Exit on error

# Define directories and files
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_FILE="${PROJECT_ROOT}/out/tables/piQTL_results_cis_trans_annotated.csv"
OUTPUT_DIR="${PROJECT_ROOT}/out/filtered_piQTL"
PYTHON_SCRIPT="${PROJECT_ROOT}/src/filter_and_plot_piQTL_by_cis.py"

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

################################################################################
# Functions
################################################################################

print_header() {
    echo ""
    echo "╔════════════════════════════════════════════════════════════════╗"
    echo "║  Filter and Plot piQTL Based on Cis-Condition                  ║"
    echo "╚════════════════════════════════════════════════════════════════╝"
    echo ""
}

print_section() {
    echo ""
    echo -e "${YELLOW}▶ $1${NC}"
    echo ""
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ $1${NC}" >&2
}

check_input_file() {
    if [[ ! -f "$INPUT_FILE" ]]; then
        print_error "Input file not found: $INPUT_FILE"
        exit 1
    fi
    print_success "Input file exists: $INPUT_FILE"
}

check_python_script() {
    if [[ ! -f "$PYTHON_SCRIPT" ]]; then
        print_error "Python script not found: $PYTHON_SCRIPT"
        exit 1
    fi
    print_success "Python script found: $PYTHON_SCRIPT"
}

check_python_dependencies() {
    print_section "Checking Python dependencies"

    python3 << EOF
import sys
try:
    import pandas
    import numpy
    import scipy
    import matplotlib
    import seaborn
    print("✓ All required Python packages are available")
except ImportError as e:
    print(f"✗ Missing package: {e}", file=sys.stderr)
    sys.exit(1)
EOF
}

run_analysis() {
    print_section "Running piQTL filtering and visualization"

    python3 "$PYTHON_SCRIPT" \
        --input "$INPUT_FILE" \
        --output "$OUTPUT_DIR"

    if [[ $? -eq 0 ]]; then
        print_success "Analysis completed successfully"
    else
        print_error "Analysis failed"
        exit 1
    fi
}

verify_outputs() {
    print_section "Verifying output files"

    local expected_files=(
        "${OUTPUT_DIR}/piQTL_results_filtered.csv"
        "${OUTPUT_DIR}/piQTL_cis_trans_with_filter.png"
        "${OUTPUT_DIR}/piQTL_cis_trans_with_filter.svg"
        "${OUTPUT_DIR}/piQTL_filtered_t_test.txt"
    )

    for file in "${expected_files[@]}"; do
        if [[ -f "$file" ]]; then
            local size=$(du -h "$file" | cut -f1)
            print_success "$(basename "$file") ($size)"
        else
            print_error "Missing output file: $file"
            exit 1
        fi
    done
}

print_summary() {
    print_section "Summary"

    echo "Output directory: $OUTPUT_DIR"
    echo ""
    echo "Generated files:"
    echo "  • piQTL_results_filtered.csv         - Filtered piQTL dataset"
    echo "  • piQTL_cis_trans_with_filter.png    - Boxplot visualization (PNG)"
    echo "  • piQTL_cis_trans_with_filter.svg    - Boxplot visualization (SVG)"
    echo "  • piQTL_filtered_t_test.txt          - Statistical test results"
    echo ""
}

################################################################################
# Main execution
################################################################################

main() {
    print_header

    print_section "Pre-flight checks"
    check_input_file
    check_python_script
    check_python_dependencies

    run_analysis
    verify_outputs

    print_summary

    echo -e "${GREEN}╔════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${GREEN}║  All tasks completed successfully!                             ║${NC}"
    echo -e "${GREEN}╚════════════════════════════════════════════════════════════════╝${NC}"
    echo ""
}

# Run main function
main
