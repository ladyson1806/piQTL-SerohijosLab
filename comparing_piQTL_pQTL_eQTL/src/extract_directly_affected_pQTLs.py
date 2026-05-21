"""Extract pQTLs which work as piQTL directory to the PPI-tagged genes."""

import argparse

import pandas as pd

# Parse command line arguments
parser = argparse.ArgumentParser()

parser.add_argument(
    "--data", type=str, help="Path to piQTLs and eQTL compared data file", dest="data"
)

parser.add_argument("--output", type=str, help="Path to output file", dest="output")

args = parser.parse_args()

# Load data
data = pd.read_csv(args.data)
# print(data.head())

# Make sure that pQTL_protein_name is str
data["pQTL_protein_name"] = data["pQTL_protein_name"].astype(str)


# parse chromosome
def is_affected(ppi: str, pqtl_genes) -> bool:
    pqtl_genes = set(pqtl_genes.split("; "))

    gene_a, gene_b = ppi.split("_")
    return gene_a in pqtl_genes or gene_b in pqtl_genes


# extract directly affected piQTLs
data["is_affected"] = data.apply(
    lambda x: is_affected(x["PPI"], x["pQTL_protein_name"]), axis=1
)

data = data[data["is_affected"]]
# print(data[["PPI", "eQTL_gene", "is_affected"]].head())
data = data.drop(columns=["is_affected"])


# save data
data.to_csv(args.output, index=False)
