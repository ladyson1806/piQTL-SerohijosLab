"""Extract overlapped pQTLs (Jakobson et al. 2025) with piQTLs."""

import argparse
from os import makedirs
from os.path import exists

import pandas as pd
from tqdm import tqdm

# Parse command line arguments
parser = argparse.ArgumentParser()

parser.add_argument("--piqtl", type=str, help="Path to piQTLs data file", dest="piqtl")

parser.add_argument("--pqtl", type=str, help="Path to pQTLs data file", dest="pqtl")

parser.add_argument(
    "--output_dir", type=str, help="Path to output directory", dest="output_dir"
)

args = parser.parse_args()

# create output directory
if not exists(args.output_dir):
    makedirs(args.output_dir)

# Load data
piqtl = pd.read_csv(args.piqtl)
pqtl = pd.read_csv(args.pqtl)
# print(piqtl.head())
# print(pqtl.head())


# merge data based on chromosome and position range
def range_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    """Check if two ranges overlap."""
    return max(a_start, b_start) <= min(a_end, b_end)


merged = []
for chr_id in tqdm(range(1, 17)):
    piqtl_chr = piqtl[piqtl["chromosome"] == chr_id]
    pqtl_chr = pqtl[pqtl["chromosome"] == chr_id]

    # remove chromosome column from eqtl_chr because it is same with piqtl_chr
    pqtl_chr = pqtl_chr.drop(columns=["chromosome"])

    for i, row_piqtl in piqtl_chr.iterrows():
        for j, row_pqtl in pqtl_chr.iterrows():
            if range_overlap(
                row_piqtl["piQTL_left"],
                row_piqtl["piQTL_right"],
                row_pqtl["pQTL_left"],
                row_pqtl["pQTL_right"],
            ):

                # check if piQTL peak is in pQTL range and vice versa
                is_piQTL_peak_in_pQTL = (
                    row_pqtl["pQTL_left"]
                    <= row_piqtl["piQTL_peak"]
                    <= row_pqtl["pQTL_right"]
                )
                is_pQTL_peak_in_piQTL = (
                    row_piqtl["piQTL_left"]
                    <= row_pqtl["pQTL_peak"]
                    <= row_piqtl["piQTL_right"]
                )

                # calculate peak distance
                peak_distance = abs(row_piqtl["piQTL_peak"] - row_pqtl["pQTL_peak"])

                # merge row from piqtl and pqtl
                merged_row = {
                    "SNP": row_piqtl["SNP"],
                    "PPI": row_piqtl["PPI"],
                    "DRUG": row_piqtl["DRUG"],
                    "pQTL_protein": row_pqtl["pQTL_protein"],
                    "pQTL_protein_name": row_pqtl["pQTL_commonName"],
                    "Chromosome": f"Chr{chr_id}",
                    "piQTL_peak": row_piqtl["piQTL_peak"],
                    "piQTL_left": row_piqtl["piQTL_left"],
                    "piQTL_right": row_piqtl["piQTL_right"],
                    "pQTL_peak": row_pqtl["pQTL_peak"],
                    "pQTL_left": row_pqtl["pQTL_left"],
                    "pQTL_right": row_pqtl["pQTL_right"],
                    "piQTL_EFFECTSIZE": row_piqtl["piQTL_EFFECTSIZE"],
                    "piQTL_pvalue": row_piqtl["piQTL_pvalue"],
                    "piQTL_FDR": row_piqtl["piQTL_FDR"],
                    "pQTL_effect": row_pqtl["pQTL_effect"],
                    "pQTL_gene_dist": row_pqtl["dist"],
                    "pQTL_cis": row_pqtl["pQTL_cis"],
                    "peak_distance": peak_distance,
                    "is_piQTL_peak_in_pQTL": is_piQTL_peak_in_pQTL,
                    "is_pQTL_peak_in_piQTL": is_pQTL_peak_in_piQTL,
                    "piQTL_locus_id": row_piqtl["piQTL_locus_id"],
                    "piQTL_locus_gene": row_piqtl["piQTL_locus_gene"],
                    "piQTL_locus_description": row_piqtl["piQTL_locus_description"],
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
    merged_df["is_piQTL_peak_in_pQTL"] & merged_df["is_pQTL_peak_in_piQTL"]
]
output_colocal = f"{args.output_dir}/piQTLs_vs_pQTLs_colocal.csv"
merged_df_colocal.to_csv(output_colocal, index=False)

# make merged dataset that the both peaks are exactly the same position
merged_df_exact = merged_df[merged_df["peak_distance"] == 0]
output_exact = f"{args.output_dir}/piQTLs_vs_pQTLs_exact.csv"
merged_df_exact.to_csv(output_exact, index=False)
