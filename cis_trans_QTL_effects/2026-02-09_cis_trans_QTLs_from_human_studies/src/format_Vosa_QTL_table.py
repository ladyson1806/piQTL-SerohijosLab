"""Formatting QTLs table from Teyssonniere's data"""

import argparse
import os

import pandas as pd


# Load the command line arguments
def load_args():
    parser = argparse.ArgumentParser(
        description="Format QTLs table from Teyssonniere's data"
    )
    parser.add_argument(
        "--input_cis",
        type=str,
        required=True,
        help="Path to the cis QTL SNP file",
    )
    parser.add_argument(
        "--input_trans",
        type=str,
        required=True,
        help="Path to the trans QTL SNP file",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to the output file",
    )
    return parser.parse_args()


# Load the QTL table
def load_QTLs(input_cis_file, input_trans_file):
    cis = pd.read_csv(input_cis_file, header=0, sep="\t")
    trans = pd.read_csv(input_trans_file, header=0, sep="\t")

    # add cis_trans column to both data
    cis["cis_trans"] = "cis"
    trans["cis_trans"] = "trans"

    # Combine cis and trans data
    QTLs = pd.concat([cis, trans], ignore_index=True)

    # Calculate absolute effect size
    QTLs["abs_z"] = QTLs["Zscore"].abs()
    return QTLs


# main function
def main():
    # Load the command line arguments
    args = load_args()

    # Load the QTLs SNPs
    QTLs = load_QTLs(args.input_cis, args.input_trans)

    # Save the QTLs to a file
    QTLs.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
