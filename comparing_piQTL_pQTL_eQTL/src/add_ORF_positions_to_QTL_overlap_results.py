"""Extract pQTLs which work as piQTL directory to the PPI-tagged genes."""

import argparse

import pandas as pd

# Parse command line arguments
parser = argparse.ArgumentParser()

parser.add_argument(
    "--overlap_input",
    type=str,
    help="Path to piQTLs, pQTLs, and eQTLs compared data file",
    dest="overlap_input",
)

parser.add_argument("--output", type=str, help="Path to output file", dest="output")

parser.add_argument(
    "--orf_positions_ref",
    type=str,
    help="Path to SGD ORF positions reference file",
    dest="orf_positions_ref",
)

parser.add_argument(
    "--target_qtl_types",
    nargs="+",
    default=[],
    help="Target QTL types to add ORF positions for (e.g., piQTL, pQTL, eQTL)",
    dest="target_qtl_types",
)

args = parser.parse_args()


# Load data
data = pd.read_csv(args.overlap_input)

# Load ORF positions reference
orf_positions = pd.read_csv(args.orf_positions_ref)
orf_positions = orf_positions[["locus_id", "name", "chrom", "start", "end"]]


# define function to add ORF positions for a given QTL type
def add_orf_positions(
    data: pd.DataFrame, ref: pd.DataFrame, qtl_type: str
) -> pd.DataFrame:
    if qtl_type == "pQTL":
        qtl_gene_col = "pQTL_protein"
    elif qtl_type == "eQTL":
        qtl_gene_col = "eQTL_gene_id"
    else:
        raise ValueError(f"Unsupported QTL type: {qtl_type}")

    # make sure the qtl_gene_col is str
    data[qtl_gene_col] = data[qtl_gene_col].astype(str)

    merged = data.merge(
        ref[["locus_id", "chrom", "start", "end"]],
        how="left",
        left_on=qtl_gene_col,
        right_on="locus_id",
        suffixes=("", f"_{qtl_type}_orf"),
    )

    merged.drop(columns=["locus_id"], inplace=True)

    merged = merged.rename(
        columns={
            "chrom": f"{qtl_type}_ORF_chrom",
            "start": f"{qtl_type}_ORF_start",
            "end": f"{qtl_type}_ORF_end",
        }
    )
    return merged


def add_orf_positions_to_piQTLs(data: pd.DataFrame, ref: pd.DataFrame) -> pd.DataFrame:
    # Since piQTLs' PPI are defined as geneA_geneB, I need to split them and add ORF positions for both genes
    data["PPI_geneA"] = data["PPI"].apply(lambda x: x.split("_")[0])
    data["PPI_geneB"] = data["PPI"].apply(lambda x: x.split("_")[1])

    # Add ORF positions for geneA
    data = data.merge(
        ref[["name", "chrom", "start", "end"]],
        how="left",
        left_on="PPI_geneA",
        right_on="name",
    )

    data.drop(columns=["PPI_geneA", "name"], inplace=True)

    data = data.rename(
        columns={
            "chrom": "piQTL_ORF_chrom_geneA",
            "start": "piQTL_ORF_start_geneA",
            "end": "piQTL_ORF_end_geneA",
        }
    )

    # Add ORF positions for geneB
    merged = data.merge(
        ref[["name", "chrom", "start", "end"]],
        how="left",
        left_on="PPI_geneB",
        right_on="name",
    )

    merged.drop(columns=["PPI_geneB", "name"], inplace=True)

    merged = merged.rename(
        columns={
            "chrom": "piQTL_ORF_chrom_geneB",
            "start": "piQTL_ORF_start_geneB",
            "end": "piQTL_ORF_end_geneB",
        }
    )

    return merged


# Add ORF positions for target QTL types
for qtl_type in args.target_qtl_types:
    if qtl_type == "piQTL":
        data = add_orf_positions_to_piQTLs(data, orf_positions)
    else:
        data = add_orf_positions(data, orf_positions, qtl_type)


# save data
data.to_csv(args.output, index=False)
