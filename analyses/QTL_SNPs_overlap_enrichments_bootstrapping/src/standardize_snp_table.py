#!/usr/bin/env python3
"""
Standardize piQTL SNP annotation table.

This script processes the piQTL_SNP_annotation.csv file and creates a standardized
table with SNP positions and their LD block boundaries.
"""

import pandas as pd
import os
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Standardize piQTL SNP annotation table."
    )
    parser.add_argument(
        "--input_file",
        type=str,
        required=True,
        help="Path to input CSV file (e.g. data/piQTL_SNP_annotation.csv)",
    )
    parser.add_argument(
        "--output_file",
        type=str,
        required=True,
        help="Path to output CSV file (e.g. out/standardized_tables/piQTL_SNP_annotation.csv)",
    )
    parser.add_argument(
        "--whitelist",
        type=str,
        default=None,
        help="Path to whitelist file containing SNP IDs to include (one SNP ID per line)",
    )
    parser.add_argument(
        "--blacklist",
        type=str,
        default=None,
        help="Path to blacklist file containing SNP IDs to exclude (one SNP ID per line)",
    )
    return parser.parse_args()


def parse_ld_block(ld_string, right_skip=0, left_skip=0):
    """
    Parse LD block string to extract first and last SNP IDs.

    Args:
        ld_string: String like "1_2" or "3" or "12_13_14_15_16_17_18_19"

    Returns:
        tuple: (first_snp_id, last_snp_id)
    """
    snps = str(ld_string).split("_")
    first_snp = int(snps[left_skip])
    last_snp = int(snps[-1 - right_skip])
    return first_snp, last_snp


def standardize_snp_table(
    input_file, output_file, whitelist_file=None, blacklist_file=None
):
    """
    Create standardized SNP table with LD block boundaries.

    Args:
        input_file: Path to input CSV file
        output_file: Path to output CSV file
    """
    # Read the input CSV
    print(f"Reading input file: {input_file}")
    df = pd.read_csv(input_file)

    print(f"Total SNPs: {len(df)}")

    # Create a dictionary mapping SNP ID to position
    snp_to_pos = dict(zip(df["SNP"], df["position"]))

    # Create the standardized table
    standardized_data = []

    for idx, row in df.iterrows():
        snp_id = row["SNP"]
        chrom = row["chrom"]
        pos = row["position"]
        ld_block = row["LD_b050"]

        # Extract only chromosome number from chrom column (e.g. "CHR_1" -> "1")
        # In addition, the mitochondrial chromosome is labeled as "MT"
        # instead of a number, so we will convert it as 17.
        chrom_int = chrom.replace("CHR_", "")
        if chrom_int == "MT":
            chrom_int = "17"

        # Get first and last SNP in LD block
        first_snp, last_snp = parse_ld_block(ld_block)

        # Check if the first and last SNPs are in the same chromosome
        is_first_same_chrom = chrom == df[df["SNP"] == first_snp]["chrom"].values[0]
        is_last_same_chrom = chrom == df[df["SNP"] == last_snp]["chrom"].values[0]
        right_skip, left_skip = 0, 0
        while not (is_first_same_chrom and is_last_same_chrom):
            # print(
            #     f"Warning: SNP {snp_id} has LD block with SNPs on different chromosomes. Adjusting LD block boundaries..."
            # )

            # Case 1: If the first SNP is not on the same chromosome,
            # skip it and check the next one in the LD block
            if not is_first_same_chrom:
                left_skip += 1
                first_snp, _ = parse_ld_block(ld_block, left_skip=left_skip)
                is_first_same_chrom = (
                    chrom == df[df["SNP"] == first_snp]["chrom"].values[0]
                )
            # Case 2: If the last SNP is not on the same chromosome,
            # skip it and check the previous one in the LD block
            elif not is_last_same_chrom:
                right_skip += 1
                _, last_snp = parse_ld_block(ld_block, right_skip=right_skip)
                is_last_same_chrom = (
                    chrom == df[df["SNP"] == last_snp]["chrom"].values[0]
                )

        # Get positions for start and end
        start = snp_to_pos[first_snp]
        end = snp_to_pos[last_snp]

        # Define SNP marker based on chromosome and position (e.g. "chr1:12345")
        snp_marker = f"chr{chrom_int}:{pos}"

        standardized_data.append(
            {
                "SNP": int(snp_id),
                "SNP_marker": snp_marker,
                "chr": int(chrom_int),
                "pos": pos,
                "start": start,
                "end": end,
            }
        )

    # Create DataFrame from standardized data
    standardized_df = pd.DataFrame(standardized_data)

    # Apply whitelist and blacklist filters if provided
    if whitelist_file:
        print(f"Applying whitelist filter: {whitelist_file}")
        with open(whitelist_file, "r") as f:
            whitelist_snps = set(int(line.strip()) for line in f if line[0] != "#")
        standardized_df = standardized_df[standardized_df["SNP"].isin(whitelist_snps)]
        print(f"SNPs after whitelist filter: {len(standardized_df)}")

    if blacklist_file:
        print(f"Applying blacklist filter: {blacklist_file}")
        with open(blacklist_file, "r") as f:
            blacklist_snps = set(int(line.strip()) for line in f if line[0] != "#")
        standardized_df = standardized_df[~standardized_df["SNP"].isin(blacklist_snps)]
        print(f"SNPs after blacklist filter: {len(standardized_df)}")

    # remove duplicates (some SNPs may appear multiple times if they are in multiple mutations)
    standardized_df = standardized_df.drop_duplicates(subset=["SNP_marker"])

    # Save to CSV
    standardized_df.to_csv(output_file, index=False)

    print(f"\n✓ Standardized table created successfully!")
    print(f"✓ Output saved to: {output_file}")
    print(f"\nFirst 10 rows of standardized table:")
    print(standardized_df.head(10).to_string(index=False))

    print(f"\nSummary statistics:")
    print(f"  - Total SNPs: {len(standardized_df)}")
    print(f"  - Chromosomes: {standardized_df['chr'].unique().tolist()}")
    print(
        f"  - LD block size range: {(standardized_df['end'] - standardized_df['start']).min()} - {(standardized_df['end'] - standardized_df['start']).max()} bp"
    )

    # check if there are negative LD block sizes (which would indicate an error in the start/end calculation)
    if (standardized_df["end"] - standardized_df["start"] < 0).any():
        print(
            f"Warning: There are negative LD block sizes, which may indicate an error in the start/end calculation."
        )
        # print the rows with negative LD block sizes for debugging
        print("\nRows with negative LD block sizes:")
        print(
            standardized_df[
                standardized_df["end"] - standardized_df["start"] < 0
            ].to_string(index=False)
        )
    else:
        print(f"✓ All LD block sizes are non-negative.")

    return standardized_df


if __name__ == "__main__":
    # Set file paths from command line arguments or use defaults
    args = parse_arguments()
    input_file = args.input_file
    output_file = args.output_file
    whitelist_file = args.whitelist
    blacklist_file = args.blacklist

    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file not found: {input_file}")
        exit(1)

    # Run standardization
    standardized_df = standardize_snp_table(
        input_file, output_file, whitelist_file, blacklist_file
    )
