"""Identifying cis pQTLs from Chick's data"""

import argparse
import os

import pandas as pd


# Load the comand line arguments
def load_args():
    parser = argparse.ArgumentParser(description="Identify cis pQTLs")
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to the pQTLs SNP file",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to the output file",
    )
    return parser.parse_args()


# Load the pQTLs SNPs
def load_pQTLs(input_file):
    pQTLs = pd.read_csv(input_file, header=0)

    # extract only significant pQTLs
    pQTLs = pQTLs[pQTLs["significant_pQTL"] == 1]

    # extract the columns of interest
    pQTLs = pQTLs[
        [
            "gene_symbol",
            "QTL_type",
            "percent_variance",
        ]
    ]

    # Calculate the beta from percent_variance
    # percent_variance = 0.5 * beta^2 * 100
    # beta = sqrt(percent_variance / 50)
    pQTLs["beta"] = (pQTLs["percent_variance"] / 50) ** 0.5

    # Add cis_trans column
    pQTLs["cis_trans"] = pQTLs["QTL_type"].apply(
        lambda x: "cis" if x == "LOCAL" else "trans"
    )
    return pQTLs


# main function
def main():
    # Load the command line arguments
    args = load_args()

    # Load the pQTLs SNPs
    pQTLs = load_pQTLs(args.input)

    # Save the pQTLs to a file
    pQTLs.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
