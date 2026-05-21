"""Add QTL information to master SNP table by merging piQTL, pQTL, and eQTL data."""

import argparse
import os
from os.path import basename, dirname, join

import pandas as pd

# Parse command line arguments
parser = argparse.ArgumentParser()

parser.add_argument(
    "--master_table",
    type=str,
    help="Path to master SNP table file",
    dest="master_snp_table",
)

parser.add_argument(
    "--piqtl", type=str, help="Path to genotype data file", dest="piqtl"
)

parser.add_argument(
    "--eqtl", type=str, help="Path to eQTLs genotype data file", dest="eqtl"
)

parser.add_argument(
    "--pqtl", type=str, help="Path to pQTLs genotype data file", dest="pqtl"
)

parser.add_argument("--output", type=str, help="Path to output file", dest="output")

args = parser.parse_args()


# Load master SNP table
master_df = pd.read_csv(args.master_snp_table)


#############################################
# Add piQTL information to master SNP table
#############################################

# Load piQTL data
piqtl = pd.read_csv(args.piqtl)

# Prepare SNP marker column for piQTL
piqtl.loc[:, "snp_marker"] = (
    "chr" + piqtl["chromosome"].astype(str) + ":" + piqtl["piQTL_peak"].astype(str)
)

# Remove unnecessary columns from piQTL
piqtl = piqtl[["snp_marker", "PPI", "DRUG"]]

# Add column indicaing piQTL target
piqtl.loc[:, "Besse_piQTL"] = piqtl["PPI"] + "(" + piqtl["DRUG"] + ")"

# Group by snp_marker to aggregate multiple targets
piqtl = (
    piqtl.groupby("snp_marker")["Besse_piQTL"]
    .apply(lambda x: ";".join(x.unique()))
    .reset_index()
)

# Merge piQTL info into master SNP table
master_df = master_df.merge(piqtl, on="snp_marker", how="left")


#############################################
# Add eQTL information to master SNP table
#############################################

# Load eQTL data
eqtl = pd.read_csv(args.eqtl)

# Prepare SNP marker column for eQTL
eqtl.loc[:, "snp_marker"] = (
    "chr" + eqtl["chromosome"].astype(str) + ":" + eqtl["eQTL_peak"].astype(str)
)

# Remove unnecessary columns from eQTL
eqtl = eqtl[["snp_marker", "eQTL_gene"]]

# Group by snp_marker to aggregate multiple targets
## Make sure the eQTL_gene column is string type
eqtl["eQTL_gene"] = eqtl["eQTL_gene"].astype(str)

eqtl = (
    eqtl.groupby("snp_marker")["eQTL_gene"]
    .apply(lambda x: ";".join(x.unique()))
    .reset_index()
)

# Merge eQTL info into master SNP table
master_df = master_df.merge(eqtl, on="snp_marker", how="left")


#############################################
# Add pQTL information to master SNP table
#############################################

# Load pQTL data
pqtl = pd.read_csv(args.pqtl)

# Prepare SNP marker column for pQTL
pqtl.loc[:, "snp_marker"] = (
    "chr" + pqtl["chromosome"].astype(str) + ":" + pqtl["pQTL_peak"].astype(str)
)

# Fill pQTL_commonName if value is empty. In that case, use pQTL_protein.
pqtl["pQTL_commonName"] = pqtl.apply(
    lambda row: (
        row["pQTL_protein"]
        if pd.isna(row["pQTL_commonName"]) or row["pQTL_commonName"] == ""
        else row["pQTL_commonName"]
    ),
    axis=1,
)

# Remove unnecessary columns from pQTL
pqtl = pqtl[["snp_marker", "pQTL_commonName"]]

# Group by snp_marker to aggregate multiple targets
pqtl = (
    pqtl.groupby("snp_marker")["pQTL_commonName"]
    .apply(lambda x: ";".join(x.unique()))
    .reset_index()
)

# Merge pQTL info into master SNP table
master_df = master_df.merge(pqtl, on="snp_marker", how="left")

# Rename columns for clarity
master_df = master_df.rename(
    columns={
        "pQTL_commonName": "Jakobson_pQTL",
        "eQTL_gene": "Albert_eQTL",
    }
)

#############################################
# Save the updated master SNP table
#############################################

# Save the updated master SNP table to output file
master_df.to_csv(args.output, index=False)

## Filter for extracting SNPs that are tested in all QTL analyses (Besse, Jakobson, Albert)
## Besse and Jakobson used indentical genotype data, but Albert used partially shared genotype samples
master_df_commonly_tested = master_df[
    ((master_df["Besse"] & master_df["Jakobson"]) & master_df["Albert"])
]
master_df_commonly_tested.to_csv(
    os.path.join(
        dirname(args.output),
        basename(args.output).replace(".csv", "_commonly_tested.csv"),
    ),
    index=False,
)


###########################################################################
# Filter for extracting SNPs that are at least one of piQTL, pQTL, or eQTL
# In other words, remove SNPs that are not piQTL, pQTL, or eQTL
###########################################################################

master_df_filtered = master_df.dropna(
    subset=["Besse_piQTL", "Jakobson_pQTL", "Albert_eQTL"], how="all"
)
master_df_filtered.to_csv(
    os.path.join(
        dirname(args.output),
        basename(args.output).replace(".csv", "_onlyQTL.csv"),
    ),
    index=False,
)

## Remove SNPs that are not piQTL, pQTL, or eQTL from the commonly tested table
master_df_commonly_tested_filtered = master_df_commonly_tested.dropna(
    subset=["Besse_piQTL", "Jakobson_pQTL", "Albert_eQTL"], how="all"
)
master_df_commonly_tested_filtered.to_csv(
    os.path.join(
        dirname(args.output),
        basename(args.output).replace(".csv", "_onlyQTL_commonly_tested.csv"),
    ),
    index=False,
)


###################################################
# Filter for extracting both of piQTL AND pQTL
###################################################

master_df_piqtl_pqtl = master_df.dropna(
    subset=["Besse_piQTL", "Jakobson_pQTL"], how="any"
)
master_df_piqtl_pqtl.to_csv(
    os.path.join(
        dirname(args.output),
        basename(args.output).replace(".csv", "_piQTL_pQTL.csv"),
    ),
    index=False,
)

## Same filter for the shared tested SNPs
master_df_commonly_tested_piqtl_pqtl = master_df_commonly_tested.dropna(
    subset=["Besse_piQTL", "Jakobson_pQTL"], how="any"
)
master_df_commonly_tested_piqtl_pqtl.to_csv(
    os.path.join(
        dirname(args.output),
        basename(args.output).replace(".csv", "_piQTL_pQTL_commonly_tested.csv"),
    ),
    index=False,
)


###################################################
# Filter for extracting both of piQTL AND eQTL
###################################################

master_df_piqtl_eqtl = master_df.dropna(
    subset=["Besse_piQTL", "Albert_eQTL"], how="any"
)
master_df_piqtl_eqtl.to_csv(
    os.path.join(
        dirname(args.output),
        basename(args.output).replace(".csv", "_piQTL_eQTL.csv"),
    ),
    index=False,
)

## Same filter for the shared tested SNPs
master_df_commonly_tested_piqtl_eqtl = master_df_commonly_tested.dropna(
    subset=["Besse_piQTL", "Albert_eQTL"], how="any"
)
master_df_commonly_tested_piqtl_eqtl.to_csv(
    os.path.join(
        dirname(args.output),
        basename(args.output).replace(".csv", "_piQTL_eQTL_commonly_tested.csv"),
    ),
    index=False,
)


###################################################
# Filter for extracting both of pQTL AND eQTL
###################################################

master_df_pqtl_eqtl = master_df.dropna(
    subset=["Jakobson_pQTL", "Albert_eQTL"], how="any"
)
master_df_pqtl_eqtl.to_csv(
    os.path.join(
        dirname(args.output),
        basename(args.output).replace(".csv", "_pQTL_eQTL.csv"),
    ),
    index=False,
)

## Same filter for the shared tested SNPs
master_df_commonly_tested_pqtl_eqtl = master_df_commonly_tested.dropna(
    subset=["Jakobson_pQTL", "Albert_eQTL"], how="any"
)
master_df_commonly_tested_pqtl_eqtl.to_csv(
    os.path.join(
        dirname(args.output),
        basename(args.output).replace(".csv", "_pQTL_eQTL_commonly_tested.csv"),
    ),
    index=False,
)


###################################################
# Filter for extracting piQTL AND pQTL AND eQTL (common to all)
###################################################

master_df_common = master_df.dropna(
    subset=["Besse_piQTL", "Jakobson_pQTL", "Albert_eQTL"], how="any"
)
master_df_common.to_csv(
    os.path.join(
        dirname(args.output),
        basename(args.output).replace(".csv", "_common_piQTL_pQTL_eQTL.csv"),
    ),
    index=False,
)

## Same filter for the shared tested SNPs
master_df_commonly_tested_common = master_df_commonly_tested.dropna(
    subset=["Besse_piQTL", "Jakobson_pQTL", "Albert_eQTL"], how="any"
)
master_df_commonly_tested_common.to_csv(
    os.path.join(
        dirname(args.output),
        basename(args.output).replace(
            ".csv", "_common_piQTL_pQTL_eQTL_commonly_tested.csv"
        ),
    ),
    index=False,
)


###################################################
# # Save summary statistics
###################################################

# Describe the common and unique SNPs between piQTL, pQTL, and eQTL
summary = {}
summary["Total_snps"] = master_df.shape[0]
summary["Total_snps_Besse_tested"] = master_df[master_df["Besse"]].shape[0]
summary["Total_snps_Jakobson_tested"] = master_df[master_df["Jakobson"]].shape[0]
summary["Total_snps_Albert_tested"] = master_df[master_df["Albert"]].shape[0]
summary["piQTL_snps"] = master_df[~master_df["Besse_piQTL"].isna()].shape[0]
summary["pQTL_snps"] = master_df[~master_df["Jakobson_pQTL"].isna()].shape[0]
summary["eQTL_snps"] = master_df[~master_df["Albert_eQTL"].isna()].shape[0]
summary["piQTL_pQTL_snps"] = master_df[
    ~master_df["Besse_piQTL"].isna() & ~master_df["Jakobson_pQTL"].isna()
].shape[0]
summary["piQTL_eQTL_snps"] = master_df[
    ~master_df["Besse_piQTL"].isna() & ~master_df["Albert_eQTL"].isna()
].shape[0]
summary["pQTL_eQTL_snps"] = master_df[
    ~master_df["Jakobson_pQTL"].isna() & ~master_df["Albert_eQTL"].isna()
].shape[0]
summary["Common_snps"] = master_df[
    ~master_df["Besse_piQTL"].isna()
    & ~master_df["Jakobson_pQTL"].isna()
    & ~master_df["Albert_eQTL"].isna()
].shape[0]

# For the commonly tested SNPs
summary["Total_snps_commonly_tested"] = master_df_commonly_tested.shape[0]
summary["piQTL_snps_commonly_tested"] = master_df_commonly_tested[
    ~master_df_commonly_tested["Besse_piQTL"].isna()
].shape[0]
summary["pQTL_snps_commonly_tested"] = master_df_commonly_tested[
    ~master_df_commonly_tested["Jakobson_pQTL"].isna()
].shape[0]
summary["eQTL_snps_commonly_tested"] = master_df_commonly_tested[
    ~master_df_commonly_tested["Albert_eQTL"].isna()
].shape[0]
summary["piQTL_pQTL_snps_commonly_tested"] = master_df_commonly_tested[
    ~master_df_commonly_tested["Besse_piQTL"].isna()
    & ~master_df_commonly_tested["Jakobson_pQTL"].isna()
].shape[0]
summary["piQTL_eQTL_snps_commonly_tested"] = master_df_commonly_tested[
    ~master_df_commonly_tested["Besse_piQTL"].isna()
    & ~master_df_commonly_tested["Albert_eQTL"].isna()
].shape[0]
summary["pQTL_eQTL_snps_commonly_tested"] = master_df_commonly_tested[
    ~master_df_commonly_tested["Jakobson_pQTL"].isna()
    & ~master_df_commonly_tested["Albert_eQTL"].isna()
].shape[0]
summary["Common_snps_commonly_tested"] = master_df_commonly_tested[
    ~master_df_commonly_tested["Besse_piQTL"].isna()
    & ~master_df_commonly_tested["Jakobson_pQTL"].isna()
    & ~master_df_commonly_tested["Albert_eQTL"].isna()
].shape[0]

# Save the summary to a text file
summary_df = pd.DataFrame(list(summary.items()), columns=["Metric", "Count"])
summary_output_path = os.path.join(dirname(args.output), "summary_QTL_snps.txt")
summary_df.to_csv(summary_output_path, index=False, sep="\t")
