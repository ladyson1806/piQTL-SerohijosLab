"""Construct combined pQTLs and eQTLs which are colocaled with piQTL."""

import argparse
import os
from os.path import basename, dirname, join

import pandas as pd

# Parse command line arguments
parser = argparse.ArgumentParser()

parser.add_argument(
    "--summary_table",
    type=str,
    help="Path to QTL comparison summary table file",
    dest="summary_table",
)

parser.add_argument(
    "--snp_annotation",
    type=str,
    help="Path to SNP annotation data file",
    dest="snp_annotation",
)

parser.add_argument(
    "--type",
    type=str,
    help="Type of summary table: exact or colocal",
    dest="type",
    choices=[
        "exact",
        "colocal_groupby_piQTL",
        "colocal_groupby_pQTL",
        "colocal_groupby_eQTL",
    ],
    default="exact",
)

parser.add_argument("--output", type=str, help="Path to output file", dest="output")

args = parser.parse_args()


###############################################
# Load and format SNP annotation table
###############################################

# Load colocaled eQTL data
snp_annotation = pd.read_csv(args.snp_annotation)

# Make snp_marker column for merging
snp_annotation.loc[:, "snp_marker"] = (
    "chr"
    + snp_annotation["chrom"].astype(str).str.replace("CHR_", "").replace("MT", "17")
    + ":"
    + snp_annotation["position"].astype(str)
)

###############################################
# Load summary table, and merge SNP annotation
###############################################

# Load summary table
summary_table = pd.read_csv(args.summary_table)

# Merge SNP annotation into summary table
summary_table_annotated = summary_table.merge(
    snp_annotation,
    on="snp_marker",
    how="left",
)


# Select and reorder columns for saving
columns_to_save_bases_front = [
    "snp_marker",
    "SNP",
    "chrom",
    "position",
    "REF",
    "ALT",
    "locus_id",
    "name",
    "sgd_id",
    "Besse",
    "Jakobson",
    "Albert",
]
columns_to_save_base_back = [
    "description",
    "snps_class_up",
    "genome_annotations",
    "snps_class_down",
    "LD_b075",
    "LD_b050",
    "is_in_gene",
    "gene_locus",
    "gene_name",
]

if args.type == "exact":
    columns_to_save_specific = [
        "Besse_piQTL",
        "Albert_eQTL",
        "Jakobson_pQTL",
    ]
elif args.type == "colocal_groupby_piQTL":
    columns_to_save_specific = [
        "PPI_DRUGs",
        "colocaled_eQTL_peaks",
        "colocaled_eQTL_genes",
        "colocaled_pQTL_peaks",
        "colocaled_pQTL_proteins",
    ]
elif args.type == "colocal_groupby_pQTL":
    columns_to_save_specific = [
        "pQTL_protein",
        "colocaled_piQTL_peaks",
        "colocaled_PPI_DRUGs",
        "colocaled_eQTL_peaks",
        "colocaled_eQTL_genes",
    ]
elif args.type == "colocal_groupby_eQTL":
    columns_to_save_specific = [
        "eQTL_gene",
        "colocaled_piQTL_peaks",
        "colocaled_PPI_DRUGs",
        "colocaled_pQTL_peaks",
        "colocaled_pQTL_proteins",
    ]

columns_to_save = (
    columns_to_save_bases_front + columns_to_save_specific + columns_to_save_base_back
)

summary_table_annotated = summary_table_annotated[columns_to_save]

# rename columns for clarity
summary_table_annotated = summary_table_annotated.rename(
    columns={
        "Besse": "isInBesse",
        "Albert": "isInAlbert",
        "Jakobson": "isInJakobson",
    }
)

# Remove duplicatio rows check by snp_marker and ALT columns
summary_table_annotated = summary_table_annotated.drop_duplicates(
    subset=["snp_marker", "ALT"]
)

# Save output
summary_table_annotated.to_csv(args.output, index=False)


###################
# Calculate basic stats
###################

total_snps = summary_table_annotated.shape[0]
besse_snps = summary_table_annotated["isInBesse"].sum()
shared_snps = summary_table_annotated[
    (summary_table_annotated["isInBesse"] == 1)
    & (summary_table_annotated["isInAlbert"] == 1)
    & (summary_table_annotated["isInJakobson"] == 1)
].shape[0]


# SNPs having three QTL effect
if args.type == "exact":
    snps_with_three_qtl = summary_table_annotated[
        (summary_table_annotated["Besse_piQTL"].notnull())
        & (summary_table_annotated["Albert_eQTL"].notnull())
        & (summary_table_annotated["Jakobson_pQTL"].notnull())
    ].shape[0]
elif args.type == "colocal_groupby_piQTL":
    snps_with_three_qtl = summary_table_annotated[
        (summary_table_annotated["PPI_DRUGs"].notnull())
        & (summary_table_annotated["colocaled_eQTL_peaks"].notnull())
        & (summary_table_annotated["colocaled_pQTL_peaks"].notnull())
    ].shape[0]
elif args.type == "colocal_groupby_pQTL":
    snps_with_three_qtl = summary_table_annotated[
        (summary_table_annotated["pQTL_protein"].notnull())
        & (summary_table_annotated["colocaled_piQTL_peaks"].notnull())
        & (summary_table_annotated["colocaled_eQTL_peaks"].notnull())
    ].shape[0]
elif args.type == "colocal_groupby_eQTL":
    snps_with_three_qtl = summary_table_annotated[
        (summary_table_annotated["eQTL_gene"].notnull())
        & (summary_table_annotated["colocaled_piQTL_peaks"].notnull())
        & (summary_table_annotated["colocaled_pQTL_peaks"].notnull())
    ].shape[0]


# save stats to a text file
stats_output_path = os.path.splitext(args.output)[0] + "_stats.txt"
with open(stats_output_path, "w") as f:
    f.write(f"Total SNPs: {total_snps}\n")
    f.write(f"SNPs in Besse et al.: {besse_snps}\n")
    f.write(f"SNPs shared among all three studies: {shared_snps}\n")
    f.write(f"SNPs with all three QTL effects: {snps_with_three_qtl}\n")
