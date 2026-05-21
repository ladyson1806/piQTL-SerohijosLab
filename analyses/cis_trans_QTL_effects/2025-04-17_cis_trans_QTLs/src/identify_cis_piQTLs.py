"""Identifying cis piQTLs"""

import argparse
import os

import pandas as pd


# Load the command line arguments
def load_args():
    parser = argparse.ArgumentParser(description="Identify cis piQTLs")
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to the piQTLs SNP file",
    )
    parser.add_argument(
        "--annotations",
        type=str,
        required=True,
        help="Path to the SNP annotation file",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to the output file",
    )
    return parser.parse_args()


# Load the piQTLs SNPs
def load_piQTLs(input_file):
    piQTLs = pd.read_csv(input_file, header=0)
    # Split PPI-formed genes into separate columns:
    # e.g. "gene1_gene2" -> "gene1", "gene2"
    piQTLs["PPI"] = piQTLs["PPI"].str.split("_")
    return piQTLs


# Load the SNP annotations
def load_annotations(annotations_file):
    annotations = pd.read_csv(annotations_file, header=0)
    # rename columns
    annotations = annotations.rename(
        columns={
            "SNP": "SNP",
            "LD050_ID": "LD050_ID",
            "cis_PPI_genes_5kb-up_5kb_down": "cis_PPI_genes",
            "cis": "cis",
        },
    )
    return annotations


# main function
def main():
    # Load the command line arguments
    args = load_args()

    # Load the piQTLs SNPs
    piQTLs = load_piQTLs(args.input)

    # Load the SNP annotations
    annotations = load_annotations(args.annotations)

    # Merge LD050_ID to the piQTLs dataframe
    piQTLs = pd.merge(piQTLs, annotations[["SNP", "LD050_ID"]], on="SNP", how="left")

    # Create table with LD050_ID and cid_PPI_genes
    ldblock_to_cisPPI = (
        annotations[["LD050_ID", "cis_PPI_genes"]].dropna().drop_duplicates()
    )

    # Merge cis-PPI genes to the LD block
    piQTLs = pd.merge(piQTLs, ldblock_to_cisPPI, on="LD050_ID", how="left")

    # Add a column to identify cis piQTLs: if cis_PPI_genes in the PPI column, it's a cis piQTL, otherwise trans
    piQTLs["cis_trans"] = piQTLs.apply(
        lambda row: "cis" if row["cis_PPI_genes"] in row["PPI"] else "trans", axis=1
    )

    # Add a column of the absolute value of the beta values
    piQTLs["abs_beta"] = piQTLs["beta"].abs()

    # Reomve outliers
    piQTLs = piQTLs[piQTLs["abs_beta"] < 1.9]

    # Extrant the columns of interest
    piQTLs = piQTLs[
        [
            "SNP",
            "beta",
            "abs_beta",
            "PPI",
            "DRUG",
            "LD050_ID",
            "cis_PPI_genes",
            "cis_trans",
        ]
    ]

    # Save the results to a file
    if os.path.dirname(args.output):
        os.makedirs(os.path.dirname(args.output), exist_ok=True)
    piQTLs.to_csv(args.output, sep="\t", index=False)
    print(f"Results saved to {args.output}")


if __name__ == "__main__":
    # Run the main function
    main()
