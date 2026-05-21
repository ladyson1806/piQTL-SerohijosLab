#!/usr/bin/env python3
"""
Filter piQTL data to PPIs with at least one cis-piQTL and create visualizations.

This script:
1. Loads piQTL results with cis/trans annotations
2. Filters to keep only PPIs that have at least one cis-piQTL
3. Computes absolute effect sizes
4. Creates a boxplot comparing cis vs trans piQTL effect sizes (colored by drug)
5. Performs independent t-test to compare cis and trans distributions
6. Saves filtered data and statistical results

Author: Analysis Pipeline
Date: 2026-02-19
"""

import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Filter piQTL data by cis-condition and create visualizations"
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input CSV file with piQTL results and cis/trans annotations",
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output directory for results"
    )
    return parser.parse_args()


def create_output_directory(output_dir):
    """Create output directory if it doesn't exist."""
    os.makedirs(output_dir, exist_ok=True)
    print(f"✓ Output directory ready: {output_dir}")


def load_and_process_data(input_file):
    """Load piQTL data and compute absolute effect sizes."""
    print(f"\nLoading data from: {input_file}")
    df = pd.read_csv(input_file)

    print(f"  Original records: {len(df)}")
    print(f"  Columns: {list(df.columns)}")

    # Compute absolute effect size
    df["abs_effectsize"] = df["piQTL_EFFECTSIZE"].abs()

    return df


def filter_to_ppi_with_cis(df):
    """Filter data to keep only PPIs that have at least one cis-piQTL."""
    print(f"\nFiltering data...")

    # Find PPIs with at least one cis-piQTL
    ppis_with_cis = set(df[df["cis_trans"] == "cis"]["PPI"].unique())
    print(f"  PPIs with at least one cis-piQTL: {len(ppis_with_cis)}")

    # Filter dataset
    df_filtered = df[df["PPI"].isin(ppis_with_cis)].copy()
    print(f"  Records after filtering: {len(df_filtered)}")
    print(f"  PPIs in filtered data: {df_filtered['PPI'].nunique()}")

    # Summary statistics of cis/trans split
    cis_count = len(df_filtered[df_filtered["cis_trans"] == "cis"])
    trans_count = len(df_filtered[df_filtered["cis_trans"] == "trans"])
    print(f"  Breakdown: {cis_count} cis, {trans_count} trans")

    return df_filtered


def save_filtered_data(df_filtered, output_dir):
    """Save filtered data to CSV."""
    output_file = os.path.join(output_dir, "piQTL_results_filtered.csv")
    df_filtered.to_csv(output_file, index=False)
    print(f"✓ Filtered data saved: {output_file}")
    return output_file


def create_boxplot(df_filtered, output_dir):
    """Create boxplot comparing cis vs trans effect sizes colored by drug."""
    print(f"\nCreating boxplot...")

    # Set up the figure
    fig, ax = plt.subplots(figsize=(12, 7))

    # Define drug colors
    drug_colors = {
        "noDrug": "#1f77b4",
        "Metformin": "#ff7f0e",
        "Fluconazole": "#2ca02c",
        "5-FC": "#d62728",
        "Trifluoperazine": "#9467bd",
    }

    sns.set_style("whitegrid")

    # Create boxplot with drug as hue
    sns.boxplot(
        data=df_filtered,
        x="cis_trans",
        y="abs_effectsize",
        hue="DRUG",
        palette=drug_colors,
        width=0.6,
        ax=ax,
    )

    # Overlay individual points with transparency (without hue to avoid legend conflicts)
    for drug in df_filtered["DRUG"].unique():
        drug_data = df_filtered[df_filtered["DRUG"] == drug]
        ax.scatter(
            [0.8 if x == "cis" else 1.8 for x in drug_data["cis_trans"]],
            drug_data["abs_effectsize"],
            alpha=0.3,
            s=40,
            color=drug_colors.get(drug, "gray"),
        )

    # Customize plot
    ax.set_xlabel("piQTL Classification", fontsize=12, fontweight="bold")
    ax.set_ylabel("|Effect Size|", fontsize=12, fontweight="bold")
    ax.set_title(
        "Distribution of piQTL Absolute Effect Sizes\n(PPIs with at least one cis-piQTL)",
        fontsize=14,
        fontweight="bold",
        pad=20,
    )
    ax.legend(
        title="Drug Condition", bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=10
    )

    plt.tight_layout()

    # Save as PNG
    png_file = os.path.join(output_dir, "piQTL_cis_trans_with_filter.png")
    fig.savefig(png_file, dpi=300, bbox_inches="tight")
    print(f"✓ Boxplot saved (PNG): {png_file}")

    # Save as SVG
    svg_file = os.path.join(output_dir, "piQTL_cis_trans_with_filter.svg")
    fig.savefig(svg_file, format="svg", bbox_inches="tight")
    print(f"✓ Boxplot saved (SVG): {svg_file}")

    plt.close(fig)


def perform_ttest(df_filtered, output_dir):
    """Perform t-test comparing cis and trans effect sizes."""
    print(f"\nPerforming statistical test...")

    # Extract absolute effect sizes for each group
    cis_effectsizes = df_filtered[df_filtered["cis_trans"] == "cis"][
        "abs_effectsize"
    ].values
    trans_effectsizes = df_filtered[df_filtered["cis_trans"] == "trans"][
        "abs_effectsize"
    ].values

    # Compute descriptive statistics
    cis_stats = {
        "n": len(cis_effectsizes),
        "mean": np.mean(cis_effectsizes),
        "std": np.std(cis_effectsizes),
        "median": np.median(cis_effectsizes),
        "min": np.min(cis_effectsizes),
        "max": np.max(cis_effectsizes),
    }

    trans_stats = {
        "n": len(trans_effectsizes),
        "mean": np.mean(trans_effectsizes),
        "std": np.std(trans_effectsizes),
        "median": np.median(trans_effectsizes),
        "min": np.min(trans_effectsizes),
        "max": np.max(trans_effectsizes),
    }

    # Perform independent samples t-test
    t_statistic, p_value = stats.ttest_ind(cis_effectsizes, trans_effectsizes)

    # Write results to file
    output_file = os.path.join(output_dir, "piQTL_filtered_t_test.txt")
    with open(output_file, "w") as f:
        f.write("=" * 70 + "\n")
        f.write("T-TEST RESULTS: cis vs trans piQTL Effect Sizes\n")
        f.write("=" * 70 + "\n\n")

        f.write("DESCRIPTIVE STATISTICS\n")
        f.write("-" * 70 + "\n")
        f.write("cis-piQTL:\n")
        f.write(f"  n              : {cis_stats['n']}\n")
        f.write(f"  Mean           : {cis_stats['mean']:.6f}\n")
        f.write(f"  Std Dev        : {cis_stats['std']:.6f}\n")
        f.write(f"  Median         : {cis_stats['median']:.6f}\n")
        f.write(f"  Min            : {cis_stats['min']:.6f}\n")
        f.write(f"  Max            : {cis_stats['max']:.6f}\n\n")

        f.write("trans-piQTL:\n")
        f.write(f"  n              : {trans_stats['n']}\n")
        f.write(f"  Mean           : {trans_stats['mean']:.6f}\n")
        f.write(f"  Std Dev        : {trans_stats['std']:.6f}\n")
        f.write(f"  Median         : {trans_stats['median']:.6f}\n")
        f.write(f"  Min            : {trans_stats['min']:.6f}\n")
        f.write(f"  Max            : {trans_stats['max']:.6f}\n\n")

        f.write("-" * 70 + "\n")
        f.write("INDEPENDENT SAMPLES T-TEST\n")
        f.write("-" * 70 + "\n")
        f.write(f"t-statistic    : {t_statistic:.6f}\n")
        f.write(f"p-value        : {p_value:.6e}\n")
        f.write(
            f"Significant    : {'Yes (p < 0.05)' if p_value < 0.05 else 'No (p >= 0.05)'}\n\n"
        )

        f.write("INTERPRETATION\n")
        f.write("-" * 70 + "\n")
        if p_value < 0.05:
            f.write("The distributions of absolute effect sizes are significantly\n")
            f.write("different between cis and trans piQTLs (p < 0.05).\n")
        else:
            f.write("No significant difference detected between cis and trans piQTL\n")
            f.write("effect size distributions (p >= 0.05).\n")
        f.write("=" * 70 + "\n")

    print(f"✓ T-test results saved: {output_file}")
    print(f"  t-statistic: {t_statistic:.6f}, p-value: {p_value:.6e}")

    return output_file


def main():
    """Main function."""
    args = parse_arguments()

    try:
        # Step 1: Create output directory
        create_output_directory(args.output)

        # Step 2: Load and process data
        df = load_and_process_data(args.input)

        # Step 3: Filter to PPIs with cis-piQTL
        df_filtered = filter_to_ppi_with_cis(df)

        # Step 4: Save filtered data
        save_filtered_data(df_filtered, args.output)

        # Step 5: Create boxplot
        create_boxplot(df_filtered, args.output)

        # Step 6: Perform t-test
        perform_ttest(df_filtered, args.output)

        print("\n" + "=" * 70)
        print("✓ ANALYSIS COMPLETE")
        print("=" * 70)
        print(f"\nOutput files saved to: {args.output}")
        print("  - piQTL_results_filtered.csv (filtered dataset)")
        print("  - piQTL_cis_trans_with_filter.png (boxplot visualization)")
        print("  - piQTL_cis_trans_with_filter.svg (boxplot visualization)")
        print("  - piQTL_filtered_t_test.txt (statistical test results)")

    except Exception as e:
        print(f"\n✗ ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
