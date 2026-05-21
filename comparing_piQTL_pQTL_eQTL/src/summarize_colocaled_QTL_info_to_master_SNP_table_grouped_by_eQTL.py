"""Construct combined pQTLs and piQTLs which are colocaled with eQTL."""

import argparse
import os
from os.path import basename, dirname, join

import pandas as pd

# Parse command line arguments
parser = argparse.ArgumentParser()

parser.add_argument(
    "--master_table", type=str, help="Path to genotype data file", dest="master_table"
)

parser.add_argument(
    "--pqtl_eqtl",
    type=str,
    help="Path to pQTLs colocalized with eQTL genotype data file",
    dest="pqtl",
)

parser.add_argument(
    "--piqtl_eqtl",
    type=str,
    help="Path to piQTLs colocalized with eQTL genotype data file",
    dest="piqtl",
)

parser.add_argument("--output", type=str, help="Path to output file", dest="output")

args = parser.parse_args()


###############################################
# Make sub-summary table for colocaled piQTLs
###############################################

# Load colocaled piQTL data
piqtl = pd.read_csv(args.piqtl)
# Extract necessary columns
piqtl = piqtl[
    [
        "PPI",
        "DRUG",
        "eQTL_gene_id",
        "eQTL_gene",
        "Chromosome",
        "piQTL_peak",
        "eQTL_peak",
    ]
].reset_index(drop=True)


# Add gene_name if protein_name is missing or empty
piqtl.loc[
    piqtl["eQTL_gene"].isna() | (piqtl["eQTL_gene"] == ""),
    "eQTL_gene",
] = (
    piqtl["eQTL_gene_id"]
    .loc[piqtl["eQTL_gene"].isna() | (piqtl["eQTL_gene"] == "")]
    .values
)

# Make SNP marker column: chr{Chromosome}:{piQTL_peak}
piqtl.loc[:, "snp_marker"] = (
    "chr"
    + piqtl["Chromosome"].astype(str).str.replace("Chr", "", regex=False)
    + ":"
    + piqtl["eQTL_peak"].astype(str)
)

# Make PPI-DRUG unique identifier column: PPI(DRUG)
piqtl.loc[:, "PPI_DRUG"] = piqtl["PPI"] + "(" + piqtl["DRUG"] + ")"

# Remove unnecessary columns
piqtl = piqtl[["snp_marker", "eQTL_gene", "PPI_DRUG", "piQTL_peak"]]

# Group by snp_marker and aggregate piQTL peaks and PPI-DRUGs
piqtl = (
    piqtl.groupby("snp_marker")
    .agg(
        {
            "PPI_DRUG": lambda x: ";".join(sorted(set(x))),
            "piQTL_peak": lambda x: ";".join(sorted(set(x.astype(str)))),
            "eQTL_gene": lambda x: ";".join(sorted(set(x))),
        }
    )
    .reset_index()
)


###############################################
# Merge sub-summary tables for colocaled pQTLs
###############################################

# Load colocaled eQTL data
pqtl = pd.read_csv(args.pqtl)

# Extract necessary columns
pqtl = pqtl[
    [
        "pQTL_protein",
        "pQTL_protein_name",
        "eQTL_gene_id",
        "eQTL_gene",
        "Chromosome",
        "pQTL_peak",
        "eQTL_peak",
    ]
].reset_index(drop=True)

# Add gene_name if protein_name is missing or empty
pqtl.loc[
    pqtl["pQTL_protein_name"].isna() | (pqtl["pQTL_protein_name"] == ""),
    "pQTL_protein_name",
] = (
    pqtl["pQTL_protein"]
    .loc[pqtl["pQTL_protein_name"].isna() | (pqtl["pQTL_protein_name"] == "")]
    .values
)

# Add gene_id if gene_name is missing or empty
pqtl.loc[pqtl["eQTL_gene"].isna() | (pqtl["eQTL_gene"] == ""), "eQTL_gene"] = (
    pqtl["eQTL_gene_id"]
    .loc[pqtl["eQTL_gene"].isna() | (pqtl["eQTL_gene"] == "")]
    .values
)

# Make SNP marker column: chr{Chromosome}:{piQTL_peak}
pqtl.loc[:, "snp_marker"] = (
    "chr"
    + pqtl["Chromosome"].astype(str).str.replace("Chr", "", regex=False)
    + ":"
    + pqtl["eQTL_peak"].astype(str)
)

# Remove unnecessary columns
pqtl = pqtl[["snp_marker", "eQTL_gene", "pQTL_peak", "pQTL_protein_name"]]

# Group by snp_marker and aggregate pQTL peaks and proteins
pqtl = (
    pqtl.groupby("snp_marker")
    .agg(
        {
            "pQTL_peak": lambda x: ";".join(sorted(set(x.astype(str)))),
            "pQTL_protein_name": lambda x: ";".join(sorted(set(x))),
            "eQTL_gene": lambda x: ";".join(sorted(set(x))),
        }
    )
    .reset_index()
)

# Rename columns to avoid confusion during merge
piqtl = piqtl.rename(
    columns={
        "piQTL_peak": "colocaled_piQTL_peaks",
        "PPI_DRUG": "colocaled_PPI_DRUGs",
        "eQTL_gene": "eQTL_gene_piQTL",
    }
)
pqtl = pqtl.rename(
    columns={
        "pQTL_peak": "colocaled_pQTL_peaks",
        "pQTL_protein_name": "colocaled_pQTL_proteins",
        "eQTL_gene": "eQTL_gene_pQTL",
    }
)

# Merge piQTL and pQTL sub-summary tables on snp_marker
colocaled_piqtl_pqtl = pd.merge(piqtl, pqtl, on="snp_marker", how="outer")

# Make eQTL_protein combined column for the merged table
colocaled_piqtl_pqtl.loc[:, "eQTL_gene"] = (
    colocaled_piqtl_pqtl["eQTL_gene_piQTL"].fillna("")
    + ";"
    + colocaled_piqtl_pqtl["eQTL_gene_pQTL"].fillna("")
)
colocaled_piqtl_pqtl.loc[:, "eQTL_gene"] = colocaled_piqtl_pqtl["eQTL_gene"].str.strip(
    ";"
)
# Remove unnecessary columns
colocaled_piqtl_pqtl = colocaled_piqtl_pqtl[
    [
        "snp_marker",
        "eQTL_gene",
        "colocaled_pQTL_peaks",
        "colocaled_pQTL_proteins",
        "colocaled_piQTL_peaks",
        "colocaled_PPI_DRUGs",
    ]
].reset_index(drop=True)

# Remove duplicated eQTL gene in eQTL_gene column
colocaled_piqtl_pqtl.loc[:, "eQTL_gene"] = colocaled_piqtl_pqtl["eQTL_gene"].apply(
    lambda x: ";".join(sorted(set(x.split(";"))))
)

# Save the merged colocaled eQTL and pQTL table
colocaled_piqtl_pqtl.to_csv(args.output, index=False)


###############################################
# Merge the colocaled sumamry table with the full genotype master SNP table
###############################################
master_snp_table = pd.read_csv(args.master_table)

# Drop unnecessary columns from master SNP table
master_snp_table = master_snp_table.drop(
    columns=["Besse_piQTL", "Jakobson_pQTL"], errors="ignore"
)

# Rename Albert_eQTL column to eQTL_gene for merging
master_snp_table = master_snp_table.rename(columns={"Albert_eQTL": "eQTL_gene"})

# From here, the full information of eQTL_gene are already in the master SNP table.
# So, drop the eQTL_gene column from colocaled_piqtl_pqtl to avoid duplication during merge
colocaled_piqtl_pqtl = colocaled_piqtl_pqtl.drop(columns=["eQTL_gene"], errors="ignore")

merged_full = pd.merge(
    master_snp_table,
    colocaled_piqtl_pqtl,
    left_on="snp_marker",
    right_on="snp_marker",
    how="outer",
)

# Save the full merged table
merged_full.to_csv(args.output.replace(".csv", "_full.csv"), index=False)
