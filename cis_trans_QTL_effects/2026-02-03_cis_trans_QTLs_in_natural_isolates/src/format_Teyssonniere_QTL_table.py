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
    parser.add_argument(
        "--target_genes",
        type=str,
        required=False,
        help="Path to the target genes file",
        default=None,
    )
    return parser.parse_args()


# Load the QTL table
def load_QTLs(input_file, target_genes=None):
    QTLs = pd.read_csv(input_file, header=0)

    # extract the columns of interest
    QTLs = QTLs[
        [
            "sid_index",
            "SNP",
            "Chr",
            "ChrPos",
            "PValue",
            "EffectSize",
            "Pheno",
            "type",
        ]
    ]

    # Calculate absolute effect size
    QTLs["abs_effect"] = QTLs["EffectSize"].abs()

    # For following plotting, add cis_trans column; valid values are {"cis" or "trans"}
    # In the original table, type column contain "TRANS" or "CIS"
    QTLs["cis_trans"] = QTLs["type"].str.lower()

    # If target_genes is provided, filter the QTLs to include only those genes
    if target_genes is not None:
        target_genes = pd.read_csv(target_genes, header=0)["ORF"]
        # In the Teyssonniere QTLs data, QTL gene name is in "Pheno" column as ORF name not gene symbol.
        QTLs = QTLs[QTLs["Pheno"].apply(lambda x: str(x) in target_genes.values)]

    return QTLs


# main function
def main():
    # Load the command line arguments
    args = load_args()

    # Load the QTLs SNPs
    QTLs = load_QTLs(args.input, args.target_genes)

    # Save the QTLs to a file
    QTLs.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
