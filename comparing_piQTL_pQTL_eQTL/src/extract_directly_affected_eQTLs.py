"""Extract eQTLs which work as pQTL directory to the pQTL genes."""

import argparse

import pandas as pd

# Parse command line arguments
parser = argparse.ArgumentParser()

parser.add_argument(
    "--data", type=str, help="Path to pQTLs and eQTL compared data file", dest="data"
)

parser.add_argument("--output", type=str, help="Path to output file", dest="output")

args = parser.parse_args()

# Load data
data = pd.read_csv(args.data)
# print(data.head())


# parse chromosome
def is_affected(eqtl_gene: str, pqtl_gene: str) -> bool:
    return eqtl_gene == pqtl_gene


# extract directly affected piQTLs
data["is_affected"] = data.apply(
    lambda x: is_affected(x["eQTL_gene_id"], x["pQTL_protein"]), axis=1
)

data = data[data["is_affected"]]
# print(data[["PPI", "eQTL_gene", "is_affected"]].head())
data = data.drop(columns=["is_affected"])


# save data
data.to_csv(args.output, index=False)
