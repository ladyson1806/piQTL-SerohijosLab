"""Extract piQTLs which work as eQTL directory to the PPI-tagged genes."""

import argparse

import pandas as pd

# Parse command line arguments
parser = argparse.ArgumentParser()

parser.add_argument(
    "--data", type=str, help="Path to piQTLs and eQTL compared data file", dest="data"
)

parser.add_argument("-o", type=str, help="Path to output file", dest="output")

args = parser.parse_args()

# Load data
data = pd.read_csv(args.data)
# print(data.head())


# parse chromosome
def is_affected(ppi: str, eqtl_gene) -> bool:
    gene_a, gene_b = ppi.split("_")
    return gene_a == eqtl_gene or gene_b == eqtl_gene


# extract directly affected piQTLs
data["is_affected"] = data.apply(
    lambda x: is_affected(x["PPI"], x["eQTL_gene"]), axis=1
)

data = data[data["is_affected"]]
# print(data[["PPI", "eQTL_gene", "is_affected"]].head())
data = data.drop(columns=["is_affected"])


# save data
data.to_csv(args.output, index=False)
