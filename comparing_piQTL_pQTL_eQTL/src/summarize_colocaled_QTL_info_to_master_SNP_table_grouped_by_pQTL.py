"""Construct combined pQTLs and eQTLs which are colocaled with pQTL."""

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
    "--piqtl_pqtl",
    type=str,
    help="Path to piQTLs colocalized with pQTL genotype data file",
    dest="piqtl",
)

parser.add_argument(
    "--pqtl_eqtl",
    type=str,
    help="Path to eQTLs colocalized with pQTL genotype data file",
    dest="eqtl",
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
        "pQTL_protein",
        "pQTL_protein_name",
        "Chromosome",
        "piQTL_peak",
        "pQTL_peak",
    ]
].reset_index(drop=True)


# Add gene_name if protein_name is missing or empty
piqtl.loc[
    piqtl["pQTL_protein_name"].isna() | (piqtl["pQTL_protein_name"] == ""),
    "pQTL_protein_name",
] = (
    piqtl["pQTL_protein"]
    .loc[piqtl["pQTL_protein_name"].isna() | (piqtl["pQTL_protein_name"] == "")]
    .values
)

# Make SNP marker column: chr{Chromosome}:{piQTL_peak}
piqtl.loc[:, "snp_marker"] = (
    "chr"
    + piqtl["Chromosome"].astype(str).str.replace("Chr", "", regex=False)
    + ":"
    + piqtl["pQTL_peak"].astype(str)
)

# Make PPI-DRUG unique identifier column: PPI(DRUG)
piqtl.loc[:, "PPI_DRUG"] = piqtl["PPI"] + "(" + piqtl["DRUG"] + ")"

# Remove unnecessary columns
piqtl = piqtl[["snp_marker", "pQTL_protein_name", "PPI_DRUG", "piQTL_peak"]]

# Group by snp_marker and aggregate piQTL peaks and PPI-DRUGs
piqtl = (
    piqtl.groupby("snp_marker")
    .agg(
        {
            "PPI_DRUG": lambda x: ";".join(sorted(set(x))),
            "piQTL_peak": lambda x: ";".join(sorted(set(x.astype(str)))),
            "pQTL_protein_name": lambda x: ";".join(sorted(set(x))),
        }
    )
    .reset_index()
)


###############################################
# Merge sub-summary tables for colocaled eQTLs
###############################################

# Load colocaled eQTL data
eqtl = pd.read_csv(args.eqtl)

# Extract necessary columns
eqtl = eqtl[
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
eqtl.loc[
    eqtl["pQTL_protein_name"].isna() | (eqtl["pQTL_protein_name"] == ""),
    "pQTL_protein_name",
] = (
    eqtl["pQTL_protein"]
    .loc[eqtl["pQTL_protein_name"].isna() | (eqtl["pQTL_protein_name"] == "")]
    .values
)

# Add gene_id if gene_name is missing or empty
eqtl.loc[eqtl["eQTL_gene"].isna() | (eqtl["eQTL_gene"] == ""), "eQTL_gene"] = (
    eqtl["eQTL_gene_id"]
    .loc[eqtl["eQTL_gene"].isna() | (eqtl["eQTL_gene"] == "")]
    .values
)

# Make SNP marker column: chr{Chromosome}:{piQTL_peak}
eqtl.loc[:, "snp_marker"] = (
    "chr"
    + eqtl["Chromosome"].astype(str).str.replace("Chr", "", regex=False)
    + ":"
    + eqtl["pQTL_peak"].astype(str)
)

# Remove unnecessary columns
eqtl = eqtl[["snp_marker", "pQTL_protein_name", "eQTL_peak", "eQTL_gene"]]

# Group by snp_marker and aggregate pQTL peaks and proteins
eqtl = (
    eqtl.groupby("snp_marker")
    .agg(
        {
            "eQTL_peak": lambda x: ";".join(sorted(set(x.astype(str)))),
            "eQTL_gene": lambda x: ";".join(sorted(set(x))),
            "pQTL_protein_name": lambda x: ";".join(sorted(set(x))),
        }
    )
    .reset_index()
)

# Rename columns to avoid confusion during merge
piqtl = piqtl.rename(
    columns={
        "piQTL_peak": "colocaled_piQTL_peaks",
        "PPI_DRUG": "colocaled_PPI_DRUGs",
        "pQTL_protein_name": "pQTL_protein_piQTL",
    }
)
eqtl = eqtl.rename(
    columns={
        "eQTL_peak": "colocaled_eQTL_peaks",
        "eQTL_gene": "colocaled_eQTL_genes",
        "pQTL_protein_name": "pQTL_protein_eQTL",
    }
)

# Merge piQTL and eQTL sub-summary tables on snp_marker
colocaled_piqtl_eqtl = pd.merge(piqtl, eqtl, on="snp_marker", how="outer")

# Make pQTL_protein combined column for the merged table
colocaled_piqtl_eqtl.loc[:, "pQTL_protein"] = (
    colocaled_piqtl_eqtl["pQTL_protein_piQTL"].fillna("")
    + ";"
    + colocaled_piqtl_eqtl["pQTL_protein_eQTL"].fillna("")
)
colocaled_piqtl_eqtl.loc[:, "pQTL_protein"] = colocaled_piqtl_eqtl[
    "pQTL_protein"
].str.strip(";")
# Remove unnecessary columns
colocaled_piqtl_eqtl = colocaled_piqtl_eqtl[
    [
        "snp_marker",
        "pQTL_protein",
        "colocaled_eQTL_peaks",
        "colocaled_eQTL_genes",
        "colocaled_piQTL_peaks",
        "colocaled_PPI_DRUGs",
    ]
].reset_index(drop=True)

# Remove duplicated pQTL protein in pQTL_protein column
colocaled_piqtl_eqtl.loc[:, "pQTL_protein"] = colocaled_piqtl_eqtl[
    "pQTL_protein"
].apply(lambda x: ";".join(sorted(set(x.split(";")))))


# Save the merged colocaled eQTL and pQTL table
colocaled_piqtl_eqtl.to_csv(args.output, index=False)


###############################################
# Merge the colocaled sumamry table with the full genotype master SNP table
###############################################
master_snp_table = pd.read_csv(args.master_table)

# Drop unnecessary columns from master SNP table
master_snp_table = master_snp_table.drop(
    columns=["Albert_eQTL", "Besse_piQTL"], errors="ignore"
)

# Rename Besse_piQTL column to PPI_DRUGs for merging
master_snp_table = master_snp_table.rename(columns={"Jakobson_pQTL": "pQTL_protein"})

# From here, the full information of pQTL_protein are already in the master SNP table.
# So, drop the pQTL_protein column from colocaled_piqtl_eqtl to avoid duplication during merge
colocaled_piqtl_eqtl = colocaled_piqtl_eqtl.drop(
    columns=["pQTL_protein"], errors="ignore"
)

merged_full = pd.merge(
    master_snp_table,
    colocaled_piqtl_eqtl,
    left_on="snp_marker",
    right_on="snp_marker",
    how="outer",
)

# Save the full merged table
merged_full.to_csv(args.output.replace(".csv", "_full.csv"), index=False)
