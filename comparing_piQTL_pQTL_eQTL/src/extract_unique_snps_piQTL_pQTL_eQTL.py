"""Extract pQTLs which work as piQTL directory to the PPI-tagged genes."""

import argparse
import os
from os.path import basename, dirname, join

import pandas as pd

# Parse command line arguments
parser = argparse.ArgumentParser()

parser.add_argument("--piqtl", type=str, help="Path to piQTLs data file", dest="piqtl")

parser.add_argument("--pqtl", type=str, help="Path to pQTLs data file", dest="pqtl")

parser.add_argument("--eqtl", type=str, help="Path to eQTLs data file", dest="eqtl")

parser.add_argument("--output", type=str, help="Path to output file", dest="output")

args = parser.parse_args()

# Prepare output directory
os.makedirs(dirname(args.output), exist_ok=True)


# Load piQTL data
piqtl = pd.read_csv(args.piqtl)
# Extract unique SNP ID and chromosome positions
piqtl_unique = piqtl[["SNP", "chromosome", "piQTL_peak"]].drop_duplicates()
# sort by SNP
piqtl_unique = piqtl_unique.sort_values(by=["SNP"])

# Write summary of unique piQTL SNPs
with open(args.output, "w") as f:
    f.write((f"Total piQTLs: {len(piqtl)}\n"))
    f.write((f"Unique piQTL SNPs: {len(piqtl_unique)}\n\n"))

# Save unique piQTL SNPs
uniq_piqtl_path = dirname(args.output) + "/piQTLs_unique_snps.csv"
piqtl_unique.to_csv(uniq_piqtl_path, index=False)


# Load eQTL data
eqtl = pd.read_csv(args.eqtl)
# Extract unique eQTL SNPs
eqtl_unique_snps = eqtl[["chromosome", "eQTL_peak"]].drop_duplicates()
# sort by chromosome and eQTL peak
eqtl_unique_snps = eqtl_unique_snps.sort_values(by=["chromosome", "eQTL_peak"])

# Write summary of unique eQTL SNPs
with open(args.output, "a") as f:
    f.write((f"Total eQTLs: {len(eqtl)}\n"))
    f.write((f"Unique eQTL SNPs: {len(eqtl_unique_snps)}\n\n"))

# Save unique eQTL SNPs
uniq_eqtl_path = dirname(args.output) + "/eQTLs_unique_snps.csv"
eqtl_unique_snps.to_csv(uniq_eqtl_path, index=False)


exit()

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
