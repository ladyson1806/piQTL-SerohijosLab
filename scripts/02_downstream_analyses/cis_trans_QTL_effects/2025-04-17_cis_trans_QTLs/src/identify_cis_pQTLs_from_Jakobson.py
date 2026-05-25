"""Identifying cis pQTLs from Chick's data"""

import argparse
import os

import pandas as pd


# Load the command line arguments
def load_args():
    parser = argparse.ArgumentParser(description="Identify cis pQTLs")
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to the pQTL SNP file",
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


# Load the pQTL SNPs
def load_pQTLs(input_file, target_genes=None):
    pQTLs = pd.read_csv(input_file, header=0)

    # remove control rows (index value is 0)
    pQTLs = pQTLs[pQTLs["index"] != 0]

    # extract the columns of interest
    pQTLs = pQTLs[
        [
            "protein",
            "commonName",
            "pVal",
            "beta",
            "isQtn",
            "index",
            "chr",
            "pos",
            "dist",
        ]
    ]

    # Convert the "Inf" in dist column to numeric infinity
    ## Before conversion, confirm the values in "dist" column that are not numeric
    ## Note: If there are any unexpected non-numeric values, this will raise an error
    pQTLs["dist"] = pQTLs["dist"].replace("Inf", float("inf"))
    pQTLs["dist"] = pd.to_numeric(pQTLs["dist"], errors="raise")

    # Add cis_trans column: trans is 1 kb or more away from the gene.
    # This value is in "dist" column.
    pQTLs["cis_trans"] = pQTLs["dist"].apply(
        lambda x: "cis" if abs(x) <= 1000 else "trans"
    )

    # Calculate absolute beta values
    pQTLs["abs_beta"] = pQTLs["beta"].abs()

    # If target_genes is provided, filter the pQTLs to include only those genes
    if target_genes is not None:
        target_genes = pd.read_csv(target_genes, header=0)["Gene"]
        # In the pQTLs data, pQTL gene name is in "commonName" column, but
        # some genes are overlapped on the genome.
        # In that case, the column have multiple genes separated by "; ",
        # such as "YPS1; YAP3".
        pQTLs = pQTLs[
            pQTLs["commonName"].apply(
                lambda x: any(
                    gene in target_genes.values for gene in str(x).split("; ")
                )
            )
        ]

    return pQTLs


# main function
def main():
    # Load the command line arguments
    args = load_args()

    # Load the pQTLs SNPs
    pQTLs = load_pQTLs(args.input, args.target_genes)

    # Save the pQTLs to a file
    pQTLs.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
