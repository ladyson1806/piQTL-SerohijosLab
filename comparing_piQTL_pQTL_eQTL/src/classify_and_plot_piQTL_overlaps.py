"""
Classify piQTLs by QTL overlaps and generate boxplots.

This script:
1. Loads piQTL effect sizes
2. Matches them to SNP overlap annotations (exact and colocalized)
3. Classifies piQTLs into 4 categories based on eQTL and pQTL overlaps
4. Generates boxplots comparing effect sizes across categories
"""

import argparse
import os
import sys

import pandas as pd
import plotly.express as px
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd


def load_args():
    """Load command line arguments."""
    parser = argparse.ArgumentParser(
        description="Classify piQTLs by QTL overlaps and generate boxplots"
    )
    parser.add_argument(
        "--piqtl_input",
        type=str,
        required=True,
        help="Path to piQTLs_formatted_lead.csv",
    )
    parser.add_argument(
        "--overlap_input",
        type=str,
        required=True,
        help="Path to overlap summary table",
    )
    parser.add_argument(
        "--output_table",
        type=str,
        required=True,
        help="Output path for summary files",
    )
    parser.add_argument(
        "--output_figures",
        type=str,
        required=True,
        help="Output directory for figures",
    )
    parser.add_argument(
        "--output_stats",
        type=str,
        required=True,
        help="Output path for statistical summary",
    )
    return parser.parse_args()


def load_piqtl_data(input_file):
    """Load piQTL effect sizes."""
    print(f"Loading piQTL data from {input_file}...")
    df = pd.read_csv(input_file)
    print(f"  Loaded {len(df)} piQTL records")
    return df


def prepare_piqtl_data(df_piqtl):
    """
    Prepare piQTL data for merging.

    - Create SNP marker: chr{chromosome}:{piQTL_peak}
    - Calculate absolute effect size
    - Return dataframe ready for merging
    """
    df = df_piqtl.copy()

    # Create SNP marker
    df["snp_marker"] = (
        "chr" + df["chromosome"].astype(str) + ":" + df["piQTL_peak"].astype(str)
    )

    # Calculate absolute effect size
    df["piQTL_EFFECTSIZE_abs"] = df["piQTL_EFFECTSIZE"].abs()

    # Extract required columns
    columns_needed = [
        "snp_marker",
        "PPI",
        "DRUG",
        "piQTL_EFFECTSIZE_abs",
    ]

    return df[columns_needed]


def load_overlap_data(input_file):
    """Load overlap summary data."""
    print(f"Loading overlap data from {input_file}...")
    df = pd.read_csv(input_file)

    # Rename columns for clarity and use common names for any overlap types
    df = df.rename(
        columns={
            "Besse_piQTL": "piQTL",  # This column is in exact match table
            "Albert_eQTL": "eQTL",  # This column is in exact match table
            "Jakobson_pQTL": "pQTL",  # This column is in exact match table
            "PPI_DRUGs": "piQTL",  # This column is in colocalized match table
            "colocaled_eQTL_genes": "eQTL",  # This column is in colocalized match table
            "colocaled_pQTL_proteins": "pQTL",  # This column is in colocalized match table
        }
    )

    # Drop rows with NaN values in piQTL column
    df = df.dropna(subset=["piQTL"])

    # check if there are duplicated rows based on snp_marker
    if not df["snp_marker"].is_unique:
        # check if the entire row is duplicated, if so, just drop duplicates
        df = df.drop_duplicates()

        # If still not unique, print warning and exit
        if not df["snp_marker"].is_unique:
            print(
                "Warning: Duplicate SNP markers found in overlap data after dropping duplicates."
            )
            print(df[df.duplicated(subset=["snp_marker"], keep=False)])
        else:
            print(
                "There are duplicated SNP markers in overlap data, but duplicates were dropped."
            )

    # Extract columns needed for follwing analysis
    columns_needed = [
        "snp_marker",
        "piQTL",
        "eQTL",
        "pQTL",
    ]

    print(f"  Loaded {len(df)} SNP records")
    return df[columns_needed]


def is_qtl_present(cell_value):
    """
    Check if a QTL is present in a cell.

    Returns True if cell is non-empty and not NaN/nan.
    """
    if pd.isna(cell_value):
        return False
    if isinstance(cell_value, str):
        cell_value = cell_value.strip()
        if cell_value == "" or cell_value.lower() == "nan":
            return False
    return True


def classify_piqtl_by_overlaps(row, eqtl_col, pqtl_col):
    """
    Classify a piQTL based on eQTL and pQTL overlap presence.

    Returns one of: "both_eQTL_pQTL", "eQTL_only", "pQTL_only", "piQTL_only"
    """
    has_eqtl = is_qtl_present(row[eqtl_col])
    has_pqtl = is_qtl_present(row[pqtl_col])

    if has_eqtl and has_pqtl:
        return "both_eQTL_pQTL"
    elif has_eqtl:
        return "eQTL_only"
    elif has_pqtl:
        return "pQTL_only"
    else:
        return "piQTL_only"


def classify_piqtl_dataset(df_piqtl, df_overlaps):
    """
    Classify piQTLs and merge with overlap data.

    Args:
        df_piqtl: piQTL dataframe with snp_marker column
        df_overlaps: Overlap summary dataframe

    Returns:
        (merged_df, summary_dict) - merged dataframe with classifications, and summary stats
    """
    # Merge piQTL data with overlap data on SNP marker
    merged = pd.merge(
        df_piqtl,
        df_overlaps,
        on="snp_marker",
        how="left",
    )

    print(merged.shape)

    # Count matches
    # full math: eQTL AND pQTL is not na
    # partial match: eQTL OR pQTL is not na
    #    partial_eQTL: eQTL is not na AND pQTL is na
    #    partial_pQTL: pQTL is not na AND eQTL is
    # unmatched: eQTL AND pQTL is na, i.e., only piQTL is present
    full_matched = merged[merged["eQTL"].notna() & merged["pQTL"].notna()].copy()
    partial_matched = merged[
        (merged["eQTL"].notna() & merged["pQTL"].isna())
        | (merged["pQTL"].notna() & merged["eQTL"].isna())
    ].copy()
    partial_eqtl = merged[merged["eQTL"].notna() & merged["pQTL"].isna()].copy()
    partial_pqtl = merged[merged["pQTL"].notna() & merged["eQTL"].isna()].copy()
    unmatched = merged[merged["eQTL"].isna() & merged["pQTL"].isna()]

    summary = ""
    summary += (
        f"  Total piQTLs: {len(merged)} (SNPs: {merged['snp_marker'].nunique()})\n"
    )
    summary += f"  Full matched piQTLs: {len(full_matched)} (SNPs: {full_matched['snp_marker'].nunique()})\n"
    summary += f"  Partial matched piQTLs: {len(partial_matched)} (SNPs: {partial_matched['snp_marker'].nunique()})\n"
    summary += f"    Partial eQTL only: {len(partial_eqtl)} (SNPs: {partial_eqtl['snp_marker'].nunique()})\n"
    summary += f"    Partial pQTL only: {len(partial_pqtl)} (SNPs: {partial_pqtl['snp_marker'].nunique()})\n"
    summary += f"  Unmatched to piQTLs: {len(unmatched)} (SNPs: {unmatched['snp_marker'].nunique()})\n"

    print(summary)

    # Classify only matched piQTLs
    merged["overlap_category"] = merged.apply(
        lambda row: classify_piqtl_by_overlaps(row, "eQTL", "pQTL"), axis=1
    )

    return merged, summary


def plot_boxplot(data, title, output_path):
    """
    Generate and save boxplot.

    Args:
        data: Dataframe with 'cis_trans' and 'piQTL_EFFECTSIZE_abs' columns
        title: Plot title
        output_path: Path for PNG file (SVG will be auto-generated)
    """
    print(f"\nGenerating boxplot: {title}")

    # Specify the category order
    category_order = ["both_eQTL_pQTL", "eQTL_only", "pQTL_only", "piQTL_only"]

    # Calculate sample counts
    sample_counts = data["overlap_category"].value_counts()
    sample_counts_text = "<br>".join(
        [f"{cat}: {sample_counts.get(cat, 0)}" for cat in category_order]
    )

    # Create boxplot
    fig = px.box(
        data,
        x="overlap_category",
        y="piQTL_EFFECTSIZE_abs",
        title=title,
        labels={
            "piQTL_EFFECTSIZE_abs": "Absolute Effect Size",
            "overlap_category": "QTL Type",
        },
        category_orders={"overlap_category": category_order},
    )

    # set y-axis to (0.2, 2.2)
    fig.update_yaxes(range=[0.2, 2.2])

    # Add sample counts as annotation
    fig.add_annotation(
        text=f"Sample counts:<br>{sample_counts_text}",
        xref="paper",
        yref="paper",
        showarrow=False,
        align="left",
        font=dict(size=11),
    )

    # Create output directory if needed
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Save as PNG
    fig.write_image(output_path)
    print(f"  Saved: {output_path}")

    # Save as SVG
    svg_path = output_path.replace(".png", ".svg")
    fig.write_image(svg_path)
    print(f"  Saved: {svg_path}")


def run_anova_and_tukey(data, output_stats_path):
    """Run one-way ANOVA and Tukey HSD on overlap categories."""
    categories = [
        "both_eQTL_pQTL",
        "eQTL_only",
        "pQTL_only",
        "piQTL_only",
    ]
    groups = []
    group_labels = []
    for category in categories:
        values = data.loc[
            data["overlap_category"] == category, "piQTL_EFFECTSIZE_abs"
        ].dropna()
        if len(values) == 0:
            continue
        groups.append(values.to_numpy())
        group_labels.append(category)

    if len(groups) < 2:
        report = "Not enough groups with data to run ANOVA.\n"
        with open(output_stats_path, "w") as f:
            f.write(report)
        return

    f_stat, p_value = stats.f_oneway(*groups)
    report_lines = [
        "ANOVA (one-way) on piQTL_EFFECTSIZE_abs by overlap_category",
        f"Groups tested: {', '.join(group_labels)}",
        f"F-statistic: {f_stat:.6g}",
        f"P-value: {p_value:.6g}",
        "",
        "Tukey HSD (alpha=0.05):",
    ]

    tukey = pairwise_tukeyhsd(
        endog=data["piQTL_EFFECTSIZE_abs"],
        groups=data["overlap_category"],
        alpha=0.05,
    )
    report_lines.append(str(tukey))
    report_lines.append("")

    with open(output_stats_path, "w") as f:
        f.write("\n".join(report_lines))


def main():
    """Main execution."""
    args = load_args()

    # Load data
    df_piqtl = load_piqtl_data(args.piqtl_input)
    df_overlap = load_overlap_data(args.overlap_input)

    # Prepare piQTL data
    df_piqtl = prepare_piqtl_data(df_piqtl)

    # Classify by overlaps
    df_classified, summary = classify_piqtl_dataset(df_piqtl, df_overlap)

    # Generate boxplots
    plot_boxplot(
        df_classified,
        "piQTL Effect Sizes by QTL Overlaps",
        args.output_figures,
    )

    # Save classification summaries
    summary_file = args.output_table.replace(
        "_summary_table.csv", "_QTL_SNP_numbers.txt"
    )
    with open(summary_file, "w") as f:
        f.write("Classification Summary:\n")
        f.write(summary)

    # Save ANOVA and Tukey test results
    run_anova_and_tukey(df_classified, args.output_stats)

    # Save classified datasets
    df_classified.to_csv(
        args.output_table,
        index=False,
    )


if __name__ == "__main__":
    main()
