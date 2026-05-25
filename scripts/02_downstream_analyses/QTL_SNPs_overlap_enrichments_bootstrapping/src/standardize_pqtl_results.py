#!/usr/bin/env python3
"""
Standardize pQTL results table.

This script processes the pQTL results CSV file and creates a standardized
table with SNP markers and pQTL boundaries.
"""

import pandas as pd
import os
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(description="Standardize pQTL results table.")
    parser.add_argument(
        "--input_file",
        type=str,
        required=True,
        help="Path to input CSV file (e.g. data/pQTL_results.csv)",
    )
    parser.add_argument(
        "--output_file",
        type=str,
        required=True,
        help="Path to output CSV file (e.g. out/standardized_tables/pQTL_results.csv)",
    )
    return parser.parse_args()


def standardize_pqtl_results(input_file, output_file):
    """
    Create standardized pQTL results table.

    Args:
        input_file: Path to input CSV file
        output_file: Path to output CSV file
    """
    # Read the input CSV
    print(f"Reading input file: {input_file}")
    df = pd.read_csv(input_file)

    print(f"Total pQTL results: {len(df)}")

    # Create the standardized table
    standardized_data = []

    for idx, row in df.iterrows():
        chr_val = int(row["chromosome"])
        pos = int(row["pQTL_peak"])
        left = int(row["pQTL_left"])
        right = int(row["pQTL_right"])

        # Ensure start <= end (in case of data entry errors)
        start = min(left, right)
        end = max(left, right)

        # Define SNP marker based on chromosome and position (e.g. "chr1:12345")
        snp_marker = f"chr{chr_val}:{pos}"

        standardized_data.append(
            {
                "SNP_marker": snp_marker,
                "chr": chr_val,
                "pos": pos,
                "start": start,
                "end": end,
            }
        )

    # Create DataFrame from standardized data
    standardized_df = pd.DataFrame(standardized_data)

    # Remove duplicates (some SNP markers may appear multiple times if they are associated with multiple proteins)
    standardized_df = standardized_df.drop_duplicates(subset=["SNP_marker"])

    # Sort by chromosome and position
    standardized_df = standardized_df.sort_values(by=["chr", "pos"]).reset_index(
        drop=True
    )

    # Save to CSV
    standardized_df.to_csv(output_file, index=False)

    print(f"\n✓ Standardized table created successfully!")
    print(f"✓ Output saved to: {output_file}")
    print(f"\nFirst 10 rows of standardized table:")
    print(standardized_df.head(10).to_string(index=False))

    print(f"\nSummary statistics:")
    print(f"  - Total pQTL SNP markers: {len(standardized_df)}")
    print(f"  - Chromosomes: {sorted(standardized_df['chr'].unique().tolist())}")
    print(
        f"  - pQTL region size range: {(standardized_df['end'] - standardized_df['start']).min()} - {(standardized_df['end'] - standardized_df['start']).max()} bp"
    )

    # Check if there are negative pQTL region sizes (which would indicate an error)
    if (standardized_df["end"] - standardized_df["start"] < 0).any():
        print(
            f"Warning: There are negative pQTL region sizes, which may indicate an error."
        )
        # Print the rows with negative sizes for debugging
        print("\nRows with negative pQTL region sizes:")
        print(
            standardized_df[
                standardized_df["end"] - standardized_df["start"] < 0
            ].to_string(index=False)
        )
    else:
        print(f"✓ All pQTL region sizes are non-negative.")

    return standardized_df


if __name__ == "__main__":
    # Set file paths from command line arguments or use defaults
    args = parse_arguments()
    input_file = args.input_file
    output_file = args.output_file

    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file not found: {input_file}")
        exit(1)

    # Run standardization
    standardized_df = standardize_pqtl_results(input_file, output_file)
