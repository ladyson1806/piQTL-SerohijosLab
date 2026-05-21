#!/usr/bin/env python3
"""
Standardize eQTL SNP annotation table.

This script processes the eQTL_SNP_annotation.csv file and creates a standardized
table with SNP markers, chromosome numbers, and positions.
"""

import argparse
import os

import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Standardize eQTL SNP annotation table."
    )
    parser.add_argument(
        "--input_file",
        type=str,
        required=True,
        help="Path to input file (e.g. data/eQTL_SNP_annotation.csv)",
    )
    parser.add_argument(
        "--output_file",
        type=str,
        required=True,
        help="Path to output CSV file (e.g. out/standardized_tables/eQTL_SNP_annotation.csv)",
    )
    return parser.parse_args()


def roman_to_int(roman):
    """
    Convert Roman numeral to integer.

    Args:
        roman: Roman numeral string (I, II, III, IV, V, VI, VII, VIII, IX, X, XI, XII, XIII, XIV, XV, XVI)

    Returns:
        int: Corresponding integer (1-16 for S. cerevisiae chromosomes)
    """
    roman_map = {
        "I": 1,
        "II": 2,
        "III": 3,
        "IV": 4,
        "V": 5,
        "VI": 6,
        "VII": 7,
        "VIII": 8,
        "IX": 9,
        "X": 10,
        "XI": 11,
        "XII": 12,
        "XIII": 13,
        "XIV": 14,
        "XV": 15,
        "XVI": 16,
    }
    return roman_map.get(roman, None)


def parse_snp_marker(snp_string):
    """
    Parse SNP marker string to extract chromosome, position, and create standardized marker.

    Args:
        snp_string: SNP string in format like "chrI:33040_A/G"

    Returns:
        tuple: (snp_marker, chr_int, position) or (None, None, None) if parsing fails
    """
    try:
        # Split by underscore to separate marker from alleles
        parts = snp_string.strip().split("_")
        marker_part = parts[0]  # e.g., "chrI:33040"

        # Extract chromosome and position from marker part
        # Format: chrX:position where X is Roman numeral
        chr_pos = marker_part.split(":")
        if len(chr_pos) != 2:
            return None, None, None

        chr_part = chr_pos[0]  # e.g., "chrI"
        pos_part = chr_pos[1]  # e.g., "33040"

        # Extract Roman numeral from chr part (remove "chr" prefix)
        roman = chr_part.replace("chr", "")

        # Convert Roman numeral to integer
        chr_int = roman_to_int(roman)
        if chr_int is None:
            return None, None, None

        # Convert position to integer
        position = int(pos_part)

        # Create standardized SNP marker (without alleles) using numeric chromosome
        # to match piQTL SNP marker format (e.g., chr1:33040 not chrI:33040)
        snp_marker = f"chr{chr_int}:{position}"

        return snp_marker, chr_int, position

    except Exception as e:
        print(f"Warning: Could not parse SNP string '{snp_string}': {e}")
        return None, None, None


def standardize_eqtl_snp_table(input_file, output_file):
    """
    Create standardized eQTL SNP table.

    Args:
        input_file: Path to input file (one SNP marker per line)
        output_file: Path to output CSV file
    """
    # Read the input file
    print(f"Reading input file: {input_file}")

    with open(input_file, "r") as f:
        snp_lines = [line.strip() for line in f if line.strip()]

    print(f"Total SNP entries: {len(snp_lines)}")

    # Parse SNP markers
    standardized_data = []

    for snp_string in snp_lines:
        snp_marker, chr_int, position = parse_snp_marker(snp_string)

        if snp_marker is not None:
            standardized_data.append(
                {
                    "SNP_marker": snp_marker,
                    "chr": chr_int,
                    "pos": position,
                }
            )

    # Create DataFrame
    standardized_df = pd.DataFrame(standardized_data)

    # Remove duplicates (in case same SNP appears multiple times)
    standardized_df = standardized_df.drop_duplicates(subset=["SNP_marker"])

    # Save to CSV
    standardized_df.to_csv(output_file, index=False)

    print(f"\n✓ Standardized table created successfully!")
    print(f"✓ Output saved to: {output_file}")
    print(f"\nFirst 10 rows of standardized table:")
    print(standardized_df.head(10).to_string(index=False))

    print(f"\nSummary statistics:")
    print(f"  - Total unique SNPs: {len(standardized_df)}")
    print(f"  - Chromosomes: {sorted(standardized_df['chr'].unique().tolist())}")
    print(
        f"  - Position range: {standardized_df['pos'].min()} - {standardized_df['pos'].max()} bp"
    )

    return standardized_df


if __name__ == "__main__":
    # Parse command line arguments
    args = parse_arguments()
    input_file = args.input_file
    output_file = args.output_file

    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file not found: {input_file}")
        exit(1)

    # Run standardization
    standardized_df = standardize_eqtl_snp_table(input_file, output_file)
