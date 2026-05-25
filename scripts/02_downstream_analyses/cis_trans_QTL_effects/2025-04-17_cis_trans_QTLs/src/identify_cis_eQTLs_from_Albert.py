"""Identifying cis eQTLs from Albert's data"""

import argparse
import os

import pandas as pd


# Load the command line arguments
def load_args():
    parser = argparse.ArgumentParser(description="Identify cis eQTLs")
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to the eQTLs SNP file",
    )
    parser.add_argument(
        "--target_genes",
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


# Load the eQTLs SNPs
def load_eQTLs(input_file):
    eQTLs = pd.read_csv(input_file, header=0)
    # extract the columns of interest
    eQTLs = eQTLs[
        [
            "gene",
            "var.exp",
            "cis",
        ]
    ]
    return eQTLs


# Load the target genes
def load_target_genes(target_genes_file):
    target_genes = pd.read_csv(target_genes_file, header=0)
    # Rename the columns
    target_genes = target_genes.rename(
        columns={"Gene": "gene_name", "ORF": "gene"},
    )
    # Extract the columns of interest
    target_genes = target_genes[["gene_name", "gene"]]
    return target_genes


# Extract only the target genes
def extract_target_genes(eQTLs, target_genes):
    # Filter the eQTLs to only include the target genes
    eQTLs = eQTLs[eQTLs["gene"].isin(target_genes["gene"])]
    # Add the gene name column
    eQTLs = pd.merge(
        eQTLs,
        target_genes,
        on="gene",
        how="left",
    )
    # Return the eQTLs with reordered columns
    eQTLs = eQTLs[["gene_name", "gene", "var.exp", "cis", "cis_trans"]]
    return eQTLs


def calculate_beta_from_var_exp(eQTLs):
    # Calculate the beta from var.exp
    # var = 0.5 x beta^2
    # beta = sqrt(var / 0.5)
    eQTLs["beta"] = (eQTLs["var.exp"] / 0.5) ** 0.5
    return eQTLs


# main function
def main():
    # Load the command line arguments
    args = load_args()

    # Load the eQTLs SNPs
    eQTLs = load_eQTLs(args.input)

    # Load the target genes
    target_genes = load_target_genes(args.target_genes)

    # Add cis_trans column
    eQTLs["cis_trans"] = eQTLs["cis"].apply(lambda x: "cis" if x else "trans")

    # Extract only the target genes' eQTLs
    eQTLs_targeted = extract_target_genes(eQTLs, target_genes)

    # Calculate the beta from var.exp
    eQTLs = calculate_beta_from_var_exp(eQTLs)
    eQTLs_targeted = calculate_beta_from_var_exp(eQTLs_targeted)

    # Save the results to a file
    eQTLs_targeted.to_csv(args.output, sep="\t", index=False)

    # Save the all eQTLs to a file
    all_eQTLs = eQTLs[
        [
            "gene",
            "var.exp",
            "beta",
            "cis",
            "cis_trans",
        ]
    ]
    all_eQTLs.to_csv(
        args.output.replace(".csv", "_all.csv"),
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    # Run the main function
    main()
