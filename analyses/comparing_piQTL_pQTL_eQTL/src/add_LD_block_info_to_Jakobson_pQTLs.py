"""Parse LB block ranges from Jakobson et al."""

import argparse

import pandas as pd

# Parse command line arguments
parser = argparse.ArgumentParser()

parser.add_argument(
    "--pQTL_input", type=str, help="Path to the pQTL input file", dest="pQTL_input"
)
parser.add_argument(
    "--ld_range_input",
    type=str,
    help="Path to LD range input file",
    dest="ld_range_input",
)

parser.add_argument("--output", type=str, help="Path to output file", dest="output")

args = parser.parse_args()

# Load Jakobson pQTL data
pQTL_df = pd.read_csv(args.pQTL_input, sep="\t")

# Load LD block range data
ld_range_df = pd.read_csv(args.ld_range_input, sep="\t")

# Extract only required columns for merging
ld_range_df = ld_range_df[["SNP", "left_pos", "right_pos"]]

# Merge pQTL data with LD block ranges on SNP ID
## piQTL data: index
## LD range data: SNP
merged_df = pd.merge(pQTL_df, ld_range_df, left_on="index", right_on="SNP", how="left")


# Rename column for saving
merged_df = merged_df.rename(
    columns={
        "chr": "chromosome",
        "abs_beta": "pQTL_effect",
        "pos": "pQTL_peak",
        "left_pos": "pQTL_left",
        "right_pos": "pQTL_right",
        "pval": "pQTL_pVal",
        "isQTN": "pQTL_isQTN",
        "protein": "pQTL_protein",
        "commonName": "pQTL_commonName",
    }
)

# Add boolean column indicating whether pQTL is cis or not
merged_df["pQTL_cis"] = merged_df["cis_trans"] == "cis"

# Save merged data with LD block info
merged_df.to_csv(args.output, sep=",", index=False)
