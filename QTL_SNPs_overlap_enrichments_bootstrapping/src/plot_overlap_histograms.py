"""
Create histograms of randomized QTL overlap counts with actual observed values highlighted.
Each figure panel has 3 subplots for pQTL, eQTL, and both overlaps.
"""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def create_histogram_panel(overlap_mode, randomized_csv, actual_csv, output_prefix):
    """
    Create a histogram panel for the specified overlap mode.

    Parameters:
    - overlap_mode: 'exact' or 'colocal'
    - randomized_csv: path to randomized overlap counts CSV
    - actual_csv: path to actual overlap status CSV
    - output_prefix: output file prefix (without extension)
    """

    # Read data
    random_df = pd.read_csv(randomized_csv)
    actual_df = pd.read_csv(actual_csv)

    # Extract actual values for this overlap mode
    actual_values = actual_df[actual_df["overlap_mode"] == overlap_mode].iloc[0]

    # Create figure with 3 subplots
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    # Color scheme
    color_pqtl = "#1f77b4"  # blue
    color_eqtl = "#ff7f0e"  # orange
    color_both = "#2ca02c"  # green
    actual_color = "#d62728"  # red

    bins = min(
        20,
        round(
            max(
                random_df["overlap_pqtl"].quantile(0.60)
                - random_df["overlap_pqtl"].min(),
                random_df["overlap_eqtl"].quantile(0.60)
                - random_df["overlap_eqtl"].min(),
                random_df["overlap_both"].quantile(0.60)
                - random_df["overlap_both"].min(),
            )
        ),
    )
    print(f"Using {bins} bins for histograms based on data range.")

    # Subplot a: pQTL
    axes[0].hist(
        random_df["overlap_pqtl"],
        bins=bins,
        color=color_pqtl,
        alpha=0.7,
        edgecolor="black",
    )
    axes[0].axvline(
        actual_values["with_pQTL"],
        color=actual_color,
        linestyle="--",
        linewidth=2.5,
        label=f"Actual: {actual_values['with_pQTL']}",
    )
    axes[0].set_xlabel("Number of Overlapping QTLs", fontsize=10)
    axes[0].set_ylabel("Frequency", fontsize=10)
    axes[0].set_title(
        f"{overlap_mode.capitalize()} Overlap: pQTL", fontsize=11, fontweight="bold"
    )
    axes[0].legend(fontsize=9)
    axes[0].grid(axis="y", alpha=0.3)

    # Subplot b: eQTL
    axes[1].hist(
        random_df["overlap_eqtl"],
        bins=bins,
        color=color_eqtl,
        alpha=0.7,
        edgecolor="black",
    )
    axes[1].axvline(
        actual_values["with_eQTL"],
        color=actual_color,
        linestyle="--",
        linewidth=2.5,
        label=f"Actual: {actual_values['with_eQTL']}",
    )
    axes[1].set_xlabel("Number of Overlapping QTLs", fontsize=10)
    axes[1].set_ylabel("Frequency", fontsize=10)
    axes[1].set_title(
        f"{overlap_mode.capitalize()} Overlap: eQTL", fontsize=11, fontweight="bold"
    )
    axes[1].legend(fontsize=9)
    axes[1].grid(axis="y", alpha=0.3)

    # Subplot c: Both
    axes[2].hist(
        random_df["overlap_both"],
        bins=bins,
        color=color_both,
        alpha=0.7,
        edgecolor="black",
    )
    axes[2].axvline(
        actual_values["both"],
        color=actual_color,
        linestyle="--",
        linewidth=2.5,
        label=f"Actual: {actual_values['both']}",
    )
    axes[2].set_xlabel("Number of Overlapping QTLs", fontsize=10)
    axes[2].set_ylabel("Frequency", fontsize=10)
    axes[2].set_title(
        f"{overlap_mode.capitalize()} Overlap: Both", fontsize=11, fontweight="bold"
    )
    axes[2].legend(fontsize=9)
    axes[2].grid(axis="y", alpha=0.3)

    plt.suptitle(
        f"QTL Overlap ({overlap_mode.upper()}): Randomized Simulation vs. Actual Observed Values",
        fontsize=12,
        fontweight="bold",
        y=0.98,
    )
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    # Save figures
    png_file = f"{output_prefix}.png"
    svg_file = f"{output_prefix}.svg"
    plt.savefig(png_file, dpi=300, bbox_inches="tight")
    plt.savefig(svg_file, bbox_inches="tight")
    print(
        f"Figures saved: {png_file} ({os.path.getsize(png_file) / 1024:.0f} KB) and {svg_file} ({os.path.getsize(svg_file) / 1024:.0f} KB)"
    )
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate histograms of randomized QTL overlap counts with actual observed values highlighted.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate exact overlap histogram
  python plot_overlap_histograms.py \\
    --overlap-mode exact \\
    --randomized-csv out/randomized_overlap_counts_exact/overlap_counts_results.csv \\
    --actual-csv data/actual_piQTL_SNPs_overlap_status.csv \\
    --output-prefix out/overlap_histograms_exact

  # Generate colocal overlap histogram
  python plot_overlap_histograms.py \\
    --overlap-mode colocal \\
    --randomized-csv out/randomized_overlap_counts_colocal/overlap_counts_results.csv \\
    --actual-csv data/actual_piQTL_SNPs_overlap_status.csv \\
    --output-prefix out/overlap_histograms_colocal
        """,
    )

    parser.add_argument(
        "--overlap-mode",
        required=True,
        choices=["exact", "colocal"],
        help="Type of overlap mode (exact or colocal)",
    )
    parser.add_argument(
        "--randomized-csv",
        required=True,
        help="Path to randomized overlap counts CSV file",
    )
    parser.add_argument(
        "--actual-csv", required=True, help="Path to actual overlap status CSV file"
    )
    parser.add_argument(
        "--output-prefix",
        required=True,
        help="Output file prefix (without extension; .png and .pdf will be added)",
    )

    args = parser.parse_args()

    # Validate input files exist
    if not os.path.exists(args.randomized_csv):
        parser.error(f"Randomized CSV file not found: {args.randomized_csv}")

    if not os.path.exists(args.actual_csv):
        parser.error(f"Actual CSV file not found: {args.actual_csv}")

    create_histogram_panel(
        args.overlap_mode, args.randomized_csv, args.actual_csv, args.output_prefix
    )
