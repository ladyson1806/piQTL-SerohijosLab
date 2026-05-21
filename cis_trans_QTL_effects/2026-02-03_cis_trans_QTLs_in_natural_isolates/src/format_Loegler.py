"""Formatting QTLs table from Loegler's data"""

import argparse
import os

import pandas as pd


# Load the command line arguments
def load_args():
    parser = argparse.ArgumentParser(
        description="Format QTLs table from Loegler's data"
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
        "--pheno_type",
        type=str,
        required=True,
        help="Specify pheno type QTLs (Valid options are 'Transcriptomics' or 'Proteomics')",
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
def load_QTLs(input_file, target_genes=None, pheno_type="Transcriptomics"):
    QTLs = pd.read_csv(input_file, header=0)

    # extract the columns of interest
    QTLs = QTLs[
        [
            "Pheno_Type",
            "Pheno",
            "Geno_Type",
            "Chr",
            "ChrPos",
            "PValue",
            "EffectSize",
            "Loc",
        ]
    ]

    # Extract only SNP rows
    QTLs = QTLs[QTLs["Geno_Type"] == "SNPs"]

    # Filter by pheno_type
    QTLs = QTLs[QTLs["Pheno_Type"] == pheno_type]

    # Format "Pheno" column. In Loegler QTLs data, "Pheno" column contain gene ORF names and QTL type,
    # such as "YCR086W.RNASeq" or "YCR088W.SMProt".
    # We only need the ORF names.
    QTLs["Pheno"] = QTLs["Pheno"].apply(lambda x: str(x).split(".")[0])

    # Calculate absolute effect size
    QTLs["abs_effect"] = QTLs["EffectSize"].abs()

    # For following plotting, add cis_trans column; valid values are {"cis" or "trans"}
    # In the original table, type column contain "Distant", "Local", or "NA".
    # We consider "Local" as "cis", "Distant" as "trans", and "NA" as None.
    def determine_cis_trans(row):
        if row["Loc"] == "Local":
            return "cis"
        elif row["Loc"] == "Distant":
            return "trans"
        else:
            return None

    QTLs["cis_trans"] = QTLs.apply(determine_cis_trans, axis=1)

    # If target_genes is provided, filter the QTLs to include only those genes
    if target_genes is not None:
        target_genes = pd.read_csv(target_genes, header=0)["ORF"]
        # In the Loegler QTLs data, QTL gene name is in "Pheno" column as ORF name not gene symbol.
        QTLs = QTLs[QTLs["Pheno"].apply(lambda x: str(x) in target_genes.values)]

    return QTLs


# main function
def main():
    # Load the command line arguments
    args = load_args()

    # Load the QTLs SNPs
    QTLs = load_QTLs(args.input, args.target_genes, args.pheno_type)

    # Save the QTLs to a file
    QTLs.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
