#!/usr/bin/env python3
"""
Add cis_trans classification to piQTL results based on proximity to PPI genes.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Add cis_trans classification to piQTL results based on proximity to PPI genes."
    )
    parser.add_argument(
        "--piqtl_results",
        type=Path,
        required=True,
        help="Path to piQTL results CSV file",
    )
    parser.add_argument(
        "--ppi_target_genes",
        type=Path,
        required=True,
        help="Path to PPI target genes CSV file with gene positions",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Path to output CSV file with cis_trans annotation",
    )
    parser.add_argument(
        "--cis-window",
        type=int,
        default=5000,
        help="Window size (bp) around gene for cis classification (default: 5000)",
    )
    return parser.parse_args()


def validate_inputs(target_genes_path, results_path):
    """Validate input files exist and have required columns."""
    if not target_genes_path.exists():
        raise FileNotFoundError(f"Target genes file not found: {target_genes_path}")
    if not results_path.exists():
        raise FileNotFoundError(f"Results file not found: {results_path}")

    # Validate required columns
    target_genes_df = pd.read_csv(target_genes_path)
    required_gene_cols = {"Gene", "Chrom", "Position_left", "Position_right"}
    if not required_gene_cols.issubset(set(target_genes_df.columns)):
        raise ValueError(
            f"Target genes file missing required columns. Expected: {required_gene_cols}"
        )

    results_df = pd.read_csv(results_path)
    required_result_cols = {"PPI", "chromosome", "piQTL_peak"}
    if not required_result_cols.issubset(set(results_df.columns)):
        raise ValueError(
            f"Results file missing required columns. Expected: {required_result_cols}"
        )

    return target_genes_df, results_df


def build_genes_dict(target_genes_df):
    """Build gene position lookup dictionary."""
    genes_dict = {}
    for idx, row in target_genes_df.iterrows():
        gene_name = row["Gene"]
        chrom = row["Chrom"]  # Format: 'CHR_12'
        pos_left = row["Position_left"]
        pos_right = row["Position_right"]
        genes_dict[(gene_name, chrom)] = (pos_left, pos_right)
    return genes_dict


def classify_cis_trans(row, genes_dict, cis_window):
    """
    Classify whether a piQTL peak is cis or trans relative to PPI genes.

    Cis: peak is within [gene_left - window, gene_right + window] range on same chromosome
    Trans: peak is far from both PPI genes
    Unknown: if PPI format invalid or genes not found
    """
    ppi_str = row["PPI"]
    chromosome = row["chromosome"]
    peak = row["piQTL_peak"]

    # Convert chromosome number to CHR format
    chrom_key = f"CHR_{chromosome}"

    # Parse PPI string (format: GeneA_GeneB)
    genes = ppi_str.split("_")
    if len(genes) != 2:
        return "unknown"

    gene_a, gene_b = genes[0], genes[1]

    # Check if peak is close to either gene
    for gene in [gene_a, gene_b]:
        gene_key = (gene, chrom_key)
        if gene_key in genes_dict:
            pos_left, pos_right = genes_dict[gene_key]
            # Define close range: left - window to right + window
            close_left = pos_left - cis_window
            close_right = pos_right + cis_window

            if close_left <= peak <= close_right:
                return "cis"

    return "trans"


def main():
    """Main function."""
    args = parse_arguments()

    # Validate inputs
    try:
        target_genes_df, results_df = validate_inputs(
            args.ppi_target_genes, args.piqtl_results
        )
    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    # Build gene dictionary
    genes_dict = build_genes_dict(target_genes_df)

    print("Loaded gene information:")
    print(f"Total genes: {len(genes_dict)}")
    print(f"Sample: {list(genes_dict.items())[:3]}")

    print(f"\nLoaded {len(results_df)} rows from {args.piqtl_results}")

    # Apply classification
    print(
        f"\nClassifying cis/trans for each row (window size: {args.cis_window} bp)..."
    )
    results_df["cis_trans"] = results_df.apply(
        lambda row: classify_cis_trans(row, genes_dict, args.cis_window), axis=1
    )

    # Print summary
    print("\nClassification summary:")
    print(results_df["cis_trans"].value_counts())

    # Create output directory if needed
    args.output.parent.mkdir(parents=True, exist_ok=True)

    # Save the updated dataframe
    results_df.to_csv(args.output, index=False)
    print(f"\nSaved updated results to {args.output}")

    # Show a sample of the updated data
    print("\nSample of updated data:")
    print(results_df[["SNP", "PPI", "chromosome", "piQTL_peak", "cis_trans"]].head(10))


if __name__ == "__main__":
    main()
