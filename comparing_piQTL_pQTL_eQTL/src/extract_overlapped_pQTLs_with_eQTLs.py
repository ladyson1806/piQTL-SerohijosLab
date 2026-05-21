"""Extract overlapped pQTLs (Jakobson et al. 2025) with eQTLs (Albert et al. 2018)."""

import argparse
from os import makedirs
from os.path import exists

import pandas as pd
from tqdm import tqdm

# Parse command line arguments
parser = argparse.ArgumentParser()

parser.add_argument("--eqtl", type=str, help="Path to eQTLs data file", dest="eqtl")

parser.add_argument("--pqtl", type=str, help="Path to pQTLs data file", dest="pqtl")

parser.add_argument(
    "--output_dir", type=str, help="Path to output directory", dest="output_dir"
)

args = parser.parse_args()

# create output directory
if not exists(args.output_dir):
    makedirs(args.output_dir)

# Load data
eqtl = pd.read_csv(args.eqtl)
pqtl = pd.read_csv(args.pqtl)
# print(eqtl.head())
# print(pqtl.head())


# merge data based on chromosome and position range
def range_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    """Check if two ranges overlap."""
    return max(a_start, b_start) <= min(a_end, b_end)


merged = []
for chr_id in tqdm(range(1, 17)):
    eqtl_chr = eqtl[eqtl["chromosome"] == chr_id]
    pqtl_chr = pqtl[pqtl["chromosome"] == chr_id]

    # remove chromosome column from eqtl_chr because it is same with piqtl_chr
    pqtl_chr = pqtl_chr.drop(columns=["chromosome"])

    for i, row_eqtl in eqtl_chr.iterrows():
        for j, row_pqtl in pqtl_chr.iterrows():
            if range_overlap(
                row_eqtl["eQTL_left"],
                row_eqtl["eQTL_right"],
                row_pqtl["pQTL_left"],
                row_pqtl["pQTL_right"],
            ):

                # check if piQTL peak is in pQTL range and vice versa
                is_eQTL_peak_in_pQTL = (
                    row_pqtl["pQTL_left"]
                    <= row_eqtl["eQTL_peak"]
                    <= row_pqtl["pQTL_right"]
                )
                is_pQTL_peak_in_eQTL = (
                    row_eqtl["eQTL_left"]
                    <= row_pqtl["pQTL_peak"]
                    <= row_eqtl["eQTL_right"]
                )

                # calculate peak distance
                peak_distance = abs(row_eqtl["eQTL_peak"] - row_pqtl["pQTL_peak"])

                # merge row from piqtl and pqtl
                merged_row = {
                    "SNP": row_pqtl["index"],
                    "pQTL_protein": row_pqtl["pQTL_protein"],
                    "pQTL_protein_name": row_pqtl["pQTL_commonName"],
                    "eQTL_gene_id": row_eqtl["eQTL_gene_id"],
                    "eQTL_gene": row_eqtl["eQTL_gene"],
                    "Chromosome": f"Chr{chr_id}",
                    "pQTL_peak": row_pqtl["pQTL_peak"],
                    "pQTL_left": row_pqtl["pQTL_left"],
                    "pQTL_right": row_pqtl["pQTL_right"],
                    "eQTL_peak": row_eqtl["eQTL_peak"],
                    "eQTL_left": row_eqtl["eQTL_left"],
                    "eQTL_right": row_eqtl["eQTL_right"],
                    "pQTL_effect": row_pqtl["pQTL_effect"],
                    "eQTL_LOD": row_eqtl["eQTL_LOD"],
                    "pQTL_gene_dist": row_pqtl["dist"],
                    "pQTL_cis": row_pqtl["pQTL_cis"],
                    "eQTL_gene_dist": row_eqtl["eQTL_gene_dist"],
                    "eQTL_cis": row_eqtl["eQTL_cis"],
                    "peak_distance": peak_distance,
                    "is_eQTL_peak_in_pQTL": is_eQTL_peak_in_pQTL,
                    "is_pQTL_peak_in_eQTL": is_pQTL_peak_in_eQTL,
                }

                # save merged row
                merged.append(merged_row)

# save merged data
merged_df = pd.DataFrame(merged)
# merged_df.to_csv(args.output.replace(".csv", "_broad.csv"), index=False)

# make partial merged data (at least one peak is in the other range)
# merged_df_partial = merged_df[
#     merged_df["is_piQTL_peak_in_eQTL"] | merged_df["is_eQTL_peak_in_piQTL"]
# ]
# merged_df_partial.to_csv(args.output.replace(".csv", "_partial.csv"), index=False)

# make co-local merged data (both peaks are in each other range)
merged_df_colocal = merged_df[
    merged_df["is_eQTL_peak_in_pQTL"] & merged_df["is_pQTL_peak_in_eQTL"]
]
output_colocal = f"{args.output_dir}/pQTLs_vs_eQTLs_colocal.csv"
merged_df_colocal.to_csv(output_colocal, index=False)

# make merged dataset that the both peaks are exactly the same position
merged_df_exact = merged_df[merged_df["peak_distance"] == 0]
output_exact = f"{args.output_dir}/pQTLs_vs_eQTLs_exact.csv"
merged_df_exact.to_csv(output_exact, index=False)
