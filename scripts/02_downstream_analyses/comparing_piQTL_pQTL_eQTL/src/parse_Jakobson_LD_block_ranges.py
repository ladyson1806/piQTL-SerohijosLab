"""Parse LB block ranges from Jakobson et al."""

import argparse

import pandas as pd

# Parse command line arguments
parser = argparse.ArgumentParser()

parser.add_argument(
    "--ld_table", type=str, help="Path to the LD block table file", dest="ld_table"
)
parser.add_argument("--ld_info", type=str, help="Path to ld info file", dest="info")

parser.add_argument("--output", type=str, help="Path to output file", dest="output")

args = parser.parse_args()

# Load LD block table
ld_table_dict = {
    "SNP": [],
    "LD_left": [],
    "LD_right": [],
}
with open(args.ld_table, "r") as f:
    for i, line in enumerate(f):
        # skip header
        if i == 0:
            continue

        # parse line
        # each line contents SNP ID in the same window, such as: 1,2,,,,,,,,
        # Note: some lines contain empty values (the comma number is adjusted to the largest window block)
        spn_window = line.strip().split(",")
        # filter empty values
        spn_window = [snp for snp in spn_window if snp != ""]
        # extract leftmost and rightmost SNPs
        left_snp = spn_window[0]
        right_snp = spn_window[-1]

        # store in dict
        ld_table_dict["SNP"].append(i)
        ld_table_dict["LD_left"].append(left_snp)
        ld_table_dict["LD_right"].append(right_snp)

ld_table = pd.DataFrame(ld_table_dict)
# print(ld_table.shape)
# print(ld_table.head())
# print(ld_table.tail())

# Make sure that the "index" column values are strings for matching
ld_table["SNP"] = ld_table["SNP"].astype(str)


# Load LD info file
ld_info = pd.read_csv(args.info, sep=",")
# print(ld_info.shape)
# print(ld_info.head())

# Make sure that the "index" column values are strings for matching
ld_info["index"] = ld_info["index"].astype(str)


# construct SNP to LD block range mapping
snp_to_ld_block = {
    "SNP": [],
    "chromosome": [],
    "pos": [],
    "left_pos": [],
    "right_pos": [],
}
for i, row in ld_table.iterrows():
    snp_id = row["SNP"]
    left_snp = row["LD_left"]
    right_snp = row["LD_right"]

    # get positions
    chromosome = ld_info[ld_info["index"] == snp_id]["chr"].values[0]
    peak_pos = ld_info[ld_info["index"] == snp_id]["pos"].values[0]
    left_pos = ld_info[ld_info["index"] == left_snp]["pos"].values[0]
    right_pos = ld_info[ld_info["index"] == right_snp]["pos"].values[0]

    # store in mapping
    snp_to_ld_block["SNP"].append(snp_id)
    snp_to_ld_block["chromosome"].append(chromosome)
    snp_to_ld_block["pos"].append(peak_pos)
    snp_to_ld_block["left_pos"].append(left_pos)
    snp_to_ld_block["right_pos"].append(right_pos)

snp_to_ld_block_df = pd.DataFrame(snp_to_ld_block)
# print(snp_to_ld_block_df.shape)
# print(snp_to_ld_block_df.head())
# print(snp_to_ld_block_df.tail())

# Save to output file
snp_to_ld_block_df.to_csv(args.output, sep="\t", index=False)
