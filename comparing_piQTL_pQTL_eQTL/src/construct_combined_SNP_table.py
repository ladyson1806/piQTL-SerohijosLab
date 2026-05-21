"""Construct combined SNP table from piQTL, pQTL, and eQTL data."""

import argparse
import os
from os.path import basename, dirname, join

import pandas as pd

# Parse command line arguments
parser = argparse.ArgumentParser()

parser.add_argument(
    "--piqtl", type=str, help="Path to genotype data file", dest="piqtl"
)

parser.add_argument(
    "--eqtl", type=str, help="Path to eQTLs genotype data file", dest="eqtl"
)

parser.add_argument("--output", type=str, help="Path to output file", dest="output")

args = parser.parse_args()


# define function to convert roman numeral to integer
def roman_to_int(s):
    roman_numerals = {
        "I": 1,
        "V": 5,
        "X": 10,
        "L": 50,
        "C": 100,
        "D": 500,
        "M": 1000,
    }
    total = 0
    prev_value = 0
    for char in reversed(s):
        value = roman_numerals[char]
        if value < prev_value:
            total -= value
        else:
            total += value
        prev_value = value
    return total


# Load piQTL data
piqtl = pd.read_csv(args.piqtl)

# Extract chromosome number and positions
piqtl_positions = piqtl[["chr", "pos"]]

# Prepare SNP marker column
piqtl_positions.loc[:, "snp_marker"] = (
    "chr"
    + piqtl_positions["chr"].astype(str)
    + ":"
    + piqtl_positions["pos"].astype(str)
)

# Add boolean column for indicating SNP usage in piQTL and pQTL analysis
piqtl_positions.loc[:, "Besse"] = True
piqtl_positions.loc[:, "Jakobson"] = True

# Load eQTL data
eqtl_dict = {
    "snp_marker": [],
    "chr": [],
    "pos": [],
}
with open(args.eqtl, "r") as f:
    for line in f:
        chr, pos_ref_alt = line.strip().split(":")
        # Remove "chr" prefix and convert roman numeral to integer
        chr_num = roman_to_int(chr.replace("chr", ""))
        # Extract position
        pos, _ = pos_ref_alt.split("_")

        # Prepare snp_marker
        snp_marker = f"chr{chr_num}:{pos}"
        eqtl_dict["snp_marker"].append(snp_marker)
        eqtl_dict["chr"].append(int(chr_num))
        eqtl_dict["pos"].append(int(pos))

eqtl_positions = pd.DataFrame(eqtl_dict)

# Add boolean column for indicating SNP usage in eQTL analysis
eqtl_positions["Albert"] = True


# Merge piQTL and eQTL data on snp_marker
merged_df = pd.merge(
    piqtl_positions[["snp_marker", "Besse", "Jakobson"]],
    eqtl_positions[["snp_marker", "Albert"]],
    on="snp_marker",
    how="outer",
)

# Fill NaN values in boolean columns with False
merged_df["Besse"] = merged_df["Besse"].fillna(False)
merged_df["Jakobson"] = merged_df["Jakobson"].fillna(False)
merged_df["Albert"] = merged_df["Albert"].fillna(False)

# Add chr and pos columns
merged_df[["chr", "pos"]] = (
    merged_df["snp_marker"]
    .str.replace("chr", "", regex=False)
    .str.split(":", expand=True)
)
merged_df["chr"] = merged_df["chr"].astype(int)
merged_df["pos"] = merged_df["pos"].astype(int)

# Sort by chr and posm
merged_df = merged_df.sort_values(by=["chr", "pos"]).reset_index(drop=True)

# Rearrange columns
merged_df = merged_df[["snp_marker", "chr", "pos", "Besse", "Jakobson", "Albert"]]

# Drop duplicate SNPs if any
merged_df = merged_df.drop_duplicates(subset=["snp_marker"])

# Save the merged data to output file
merged_df.to_csv(args.output, index=False)


# Describe the common and unique SNPs between piQTL&pQTL and eQTL
summary = {}
summary["Total_snps"] = merged_df.shape[0]
summary["piQTL_pQTL_snps"] = merged_df[
    (merged_df["Besse"]) | (merged_df["Jakobson"])
].shape[0]
summary["eQTL_snps"] = merged_df[merged_df["Albert"]].shape[0]
summary["Common_snps"] = merged_df[
    (merged_df["Albert"]) & ((merged_df["Besse"]) | (merged_df["Jakobson"]))
].shape[0]

# Save the summary to a text file
summary_df = pd.DataFrame(list(summary.items()), columns=["Metric", "Count"])
summary_output_path = os.path.join(dirname(args.output), "summary_common_snps.txt")
summary_df.to_csv(summary_output_path, index=False, sep="\t")
print("Summary of common and unique SNPs saved to:", summary_output_path)
