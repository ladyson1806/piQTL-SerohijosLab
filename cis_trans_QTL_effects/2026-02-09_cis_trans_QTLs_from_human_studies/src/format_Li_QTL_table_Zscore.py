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
        "--input",
        type=str,
        required=True,
        help="Path to the QTL SNP file",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to the output file",
    )
    return parser.parse_args()


# Load the QTL table
def load_QTLs(input_file):
    QTLs = pd.read_csv(input_file, header=0)

    # extract the columns of interest
    QTLs = QTLs[
        [
            "BETA",
            "cis_trans",
        ]
    ]

    # Calxulate Zscore from BETA
    QTLs["Zscore"] = (QTLs["BETA"] - QTLs["BETA"].mean()) / QTLs["BETA"].std()

    # Calculate absolute effect size
    QTLs["abs_z"] = QTLs["Zscore"].abs()
    return QTLs


# main function
def main():
    # Load the command line arguments
    args = load_args()

    # Load the QTLs SNPs
    QTLs = load_QTLs(args.input)

    # Save the QTLs to a file
    QTLs.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
