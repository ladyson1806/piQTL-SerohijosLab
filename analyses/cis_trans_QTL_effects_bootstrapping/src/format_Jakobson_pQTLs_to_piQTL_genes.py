"""Filter Jakobson pQTLs to only target piQTL genes"""

import argparse
import os

import pandas as pd


def load_args():
    """Load command line arguments"""
    parser = argparse.ArgumentParser(
        description="Filter Jakobson pQTLs to only target piQTL genes"
    )
    parser.add_argument(
        "--input_pqtl",
        type=str,
        required=True,
        help="Path to the Jakobson pQTLs input file",
    )
    parser.add_argument(
        "--input_target_genes",
        type=str,
        required=True,
        help="Path to the target piQTL genes file",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to the output file",
    )
    return parser.parse_args()


def load_target_genes(target_genes_file):
    """Load target gene ORFs from the target genes file"""
    target_genes_df = pd.read_csv(target_genes_file)

    # Extract the ORF column as a set for efficient lookup
    target_orfs = set(target_genes_df["ORF"].values)

    print(f"Loaded {len(target_orfs)} target genes")
    return target_orfs


def load_pqtls(pqtl_file):
    """Load pQTL data"""
    pqtls = pd.read_csv(pqtl_file, sep="\t")
    print(f"Loaded {len(pqtls)} pQTL records")
    return pqtls


def filter_pqtls(pqtls, target_orfs):
    """Filter pQTLs to only include target genes"""
    # Filter where protein column matches any target ORF
    filtered = pqtls[pqtls["protein"].isin(target_orfs)].copy()
    print(f"Filtered to {len(filtered)} pQTL records targeting piQTL genes")
    return filtered


def main():
    # Load arguments
    args = load_args()

    # Load target genes
    target_orfs = load_target_genes(args.input_target_genes)

    # Load pQTL data
    pqtls = load_pqtls(args.input_pqtl)

    # Filter pQTLs
    filtered_pqtls = filter_pqtls(pqtls, target_orfs)

    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    # Save filtered pQTLs
    filtered_pqtls.to_csv(args.output, sep="\t", index=False)
    print(f"Saved filtered pQTLs to {args.output}")


if __name__ == "__main__":
    main()
