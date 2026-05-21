#!/usr/bin/env python3
"""
Extract piQTL SNPs that are also present in eQTL data.

This script identifies shared SNPs between piQTL and eQTL datasets
and creates a filtered piQTL table containing only the shared SNPs.
"""

import argparse
import os

import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(description="Extract piQTL-eQTL shared SNPs.")
    parser.add_argument(
        "--piqtl_snp_annotation",
        type=str,
        required=True,
        help="Path to standardized piQTL SNP annotation file (e.g. out/standardized_tables/piQTL_SNP_annotation.csv)",
    )
    parser.add_argument(
        "--eqtl_snp_annotation",
        type=str,
        required=True,
        help="Path to standardized eQTL SNP annotation file (e.g. out/standardized_tables/eQTL_SNP_annotation.csv)",
    )
    parser.add_argument(
        "--output_file",
        type=str,
        required=True,
        help="Path to output CSV file (e.g. out/standardized_tables/piQTL_SNP_annotation_shared_with_eQTLs.csv)",
    )
    return parser.parse_args()


def extract_shared_snps(piqtl_file, eqtl_file, output_file):
    """
    Extract piQTL SNPs that are also present in eQTL data.

    Args:
        piqtl_file: Path to standardized piQTL SNP annotation CSV
        eqtl_file: Path to standardized eQTL SNP annotation CSV
        output_file: Path to output CSV file
    """
    # Read the piQTL SNP annotation file
    print(f"Reading piQTL SNP annotation file: {piqtl_file}")
    piqtl_df = pd.read_csv(piqtl_file)
    print(f"Total piQTL SNPs: {len(piqtl_df)}")

    # Read the eQTL SNP annotation file
    print(f"Reading eQTL SNP annotation file: {eqtl_file}")
    eqtl_df = pd.read_csv(eqtl_file)
    print(f"Total eQTL SNPs: {len(eqtl_df)}")

    # Extract SNP markers from eQTL data
    eqtl_snp_markers = set(eqtl_df["SNP_marker"].unique())
    print(f"Unique eQTL SNP markers: {len(eqtl_snp_markers)}")

    # Filter piQTL SNPs to keep only those present in eQTL data
    shared_snps_df = piqtl_df[piqtl_df["SNP_marker"].isin(eqtl_snp_markers)]

    print(f"\n✓ Shared SNPs found: {len(shared_snps_df)}")

    # Save to CSV
    shared_snps_df.to_csv(output_file, index=False)

    print(f"✓ Output saved to: {output_file}")
    print(f"\nFirst 10 rows of shared SNP table:")
    print(shared_snps_df.head(10).to_string(index=False))

    print(f"\nSummary statistics:")
    print(f"  - Total shared SNPs: {len(shared_snps_df)}")
    print(f"  - Chromosomes: {sorted(shared_snps_df['chr'].unique().tolist())}")
    print(
        f"  - Position range: {shared_snps_df['pos'].min()} - {shared_snps_df['pos'].max()} bp"
    )

    # Calculate coverage statistics
    coverage_pct = (len(shared_snps_df) / len(piqtl_df)) * 100
    print(
        f"  - piQTL SNPs with eQTL support: {len(shared_snps_df)} / {len(piqtl_df)} ({coverage_pct:.2f}%)"
    )

    return shared_snps_df


if __name__ == "__main__":
    # Parse command line arguments
    args = parse_arguments()
    piqtl_file = args.piqtl_snp_annotation
    eqtl_file = args.eqtl_snp_annotation
    output_file = args.output_file

    # Check if input files exist
    if not os.path.exists(piqtl_file):
        print(f"Error: piQTL SNP annotation file not found: {piqtl_file}")
        exit(1)

    if not os.path.exists(eqtl_file):
        print(f"Error: eQTL SNP annotation file not found: {eqtl_file}")
        exit(1)

    # Extract shared SNPs
    shared_snps_df = extract_shared_snps(piqtl_file, eqtl_file, output_file)
