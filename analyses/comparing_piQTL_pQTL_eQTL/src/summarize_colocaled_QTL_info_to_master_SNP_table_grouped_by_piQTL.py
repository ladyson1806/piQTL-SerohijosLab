"""Construct combined pQTLs and eQTLs which are colocaled with piQTL."""

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
    "--piqtl_eqtl",
    type=str,
    help="Path to eQTLs colocalized with piQTL genotype data file",
    dest="eqtl",
)

parser.add_argument(
    "--piqtl_pqtl",
    type=str,
    help="Path to pQTLs colocalized with pQTL genotype data file",
    dest="pqtl",
)

parser.add_argument("--output", type=str, help="Path to output file", dest="output")

args = parser.parse_args()


###############################################
# Make sub-summary table for colocaled eQTLs
###############################################

# Load colocaled eQTL data
eqtl = pd.read_csv(args.eqtl)

# Extract necessary columns
eqtl = eqtl[
    [
        "PPI",
        "DRUG",
        "eQTL_gene",
        "eQTL_gene_id",
        "Chromosome",
        "piQTL_peak",
        "eQTL_peak",
    ]
].reset_index(drop=True)

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
    + eqtl["piQTL_peak"].astype(str)
)

# Make PPI-DRUG unique identifier column: PPI(DRUG)
eqtl.loc[:, "PPI_DRUG"] = eqtl["PPI"] + "(" + eqtl["DRUG"] + ")"

# Remove unnecessary columns
eqtl = eqtl[["snp_marker", "PPI_DRUG", "eQTL_peak", "eQTL_gene"]]

# Group by snp_marker and aggregate eQTL peaks and genes
eqtl = (
    eqtl.groupby("snp_marker")
    .agg(
        {
            "PPI_DRUG": lambda x: ";".join(sorted(set(x))),
            "eQTL_peak": lambda x: ";".join(sorted(set(x.astype(str)))),
            "eQTL_gene": lambda x: ";".join(sorted(set(x))),
        }
    )
    .reset_index()
)


###############################################
# Merge sub-summary tables for colocaled pQTLs
###############################################

# Load colocaled pQTL data
pqtl = pd.read_csv(args.pqtl)

# Extract necessary columns
pqtl = pqtl[
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
pqtl.loc[
    pqtl["pQTL_protein_name"].isna() | (pqtl["pQTL_protein_name"] == ""),
    "pQTL_protein_name",
] = (
    pqtl["pQTL_protein"]
    .loc[pqtl["pQTL_protein_name"].isna() | (pqtl["pQTL_protein_name"] == "")]
    .values
)

# Make SNP marker column: chr{Chromosome}:{piQTL_peak}
pqtl.loc[:, "snp_marker"] = (
    "chr"
    + pqtl["Chromosome"].astype(str).str.replace("Chr", "", regex=False)
    + ":"
    + pqtl["piQTL_peak"].astype(str)
)

# Make PPI-DRUG unique identifier column: PPI(DRUG)
pqtl.loc[:, "PPI_DRUG"] = pqtl["PPI"] + "(" + pqtl["DRUG"] + ")"

# Remove unnecessary columns
pqtl = pqtl[["snp_marker", "PPI_DRUG", "pQTL_peak", "pQTL_protein_name"]]

# Group by snp_marker and aggregate pQTL peaks and proteins
pqtl = (
    pqtl.groupby("snp_marker")
    .agg(
        {
            "PPI_DRUG": lambda x: ";".join(sorted(set(x))),
            "pQTL_peak": lambda x: ";".join(sorted(set(x.astype(str)))),
            "pQTL_protein_name": lambda x: ";".join(sorted(set(x))),
        }
    )
    .reset_index()
)

# Rename columns to avoid confusion during merge
eqtl = eqtl.rename(
    columns={
        "eQTL_peak": "colocaled_eQTL_peaks",
        "eQTL_gene": "colocaled_eQTL_genes",
        "PPI_DRUG": "eQTL_PPI_DRUGs",
    }
)
pqtl = pqtl.rename(
    columns={
        "pQTL_peak": "colocaled_pQTL_peaks",
        "pQTL_protein_name": "colocaled_pQTL_proteins",
        "PPI_DRUG": "pQTL_PPI_DRUGs",
    }
)

# Merge eQTL and pQTL sub-summary tables on snp_marker
colocaled_eqtl_pqtl = pd.merge(eqtl, pqtl, on="snp_marker", how="outer")

# Make PPI_DRUG combined column for the merged table
colocaled_eqtl_pqtl.loc[:, "PPI_DRUGs"] = (
    colocaled_eqtl_pqtl["eQTL_PPI_DRUGs"].fillna("")
    + ";"
    + colocaled_eqtl_pqtl["pQTL_PPI_DRUGs"].fillna("")
)
colocaled_eqtl_pqtl.loc[:, "PPI_DRUGs"] = colocaled_eqtl_pqtl["PPI_DRUGs"].str.strip(
    ";"
)
# Remove unnecessary columns
colocaled_eqtl_pqtl = colocaled_eqtl_pqtl[
    [
        "snp_marker",
        "PPI_DRUGs",
        "colocaled_eQTL_peaks",
        "colocaled_eQTL_genes",
        "colocaled_pQTL_peaks",
        "colocaled_pQTL_proteins",
    ]
].reset_index(drop=True)

# Remove duplicated PPI_DRUG in PPI_DRUGs column
colocaled_eqtl_pqtl.loc[:, "PPI_DRUGs"] = colocaled_eqtl_pqtl["PPI_DRUGs"].apply(
    lambda x: ";".join(sorted(set(x.split(";"))))
)


# Save the merged colocaled eQTL and pQTL table
colocaled_eqtl_pqtl.to_csv(args.output, index=False)


###############################################
# Merge the colocaled sumamry table with the full genotype master SNP table
###############################################
master_snp_table = pd.read_csv(args.master_table)

# Drop unnecessary columns from master SNP table
master_snp_table = master_snp_table.drop(
    columns=["Albert_eQTL", "Jakobson_pQTL"], errors="ignore"
)

# Rename Besse_piQTL column to PPI_DRUGs for merging
master_snp_table = master_snp_table.rename(columns={"Besse_piQTL": "PPI_DRUGs"})

# From here, the full information of PPI_DURGs are already in the master SNP table.
# So, drop the PPI_DRUGs column from colocaled_eqtl_pqtl to avoid duplication during merge
colocaled_eqtl_pqtl = colocaled_eqtl_pqtl.drop(columns=["PPI_DRUGs"], errors="ignore")

merged_full = pd.merge(
    master_snp_table,
    colocaled_eqtl_pqtl,
    left_on="snp_marker",
    right_on="snp_marker",
    how="outer",
)

# Save the full merged table
merged_full.to_csv(args.output.replace(".csv", "_full.csv"), index=False)
