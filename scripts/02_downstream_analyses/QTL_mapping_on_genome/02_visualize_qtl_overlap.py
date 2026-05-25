#!/usr/bin/env python3
"""
Visualize QTL overlap patterns across the genome.

Creates a multi-panel figure with:
1. Genome-wide heatmap of QTL counts
2. Albert SNP indicator strip
3-5. Three scatter plots for piQTL, pQTL, and eQTL counts

Usage:
    python visualize_qtl_overlap.py --mode exact
    python visualize_qtl_overlap.py --mode colocal --output-dir figures/
"""

import argparse
import warnings
from pathlib import Path

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm
from matplotlib.ticker import ScalarFormatter

warnings.filterwarnings("ignore")

# Constants
CM = 1 / 2.54  # Convert cm to inches for matplotlib


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Visualize QTL overlap patterns")
    parser.add_argument(
        "--mode",
        type=str,
        required=True,
        choices=["exact", "colocal"],
        help="QTL overlap mode: exact (peak=position) or colocal (position in range)",
    )
    parser.add_argument(
        "--input",
        type=str,
        default="out/tables/QTL_overlap_summary.csv",
        help="Input CSV file with QTL overlap summary",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="out/figures",
        help="Output directory for figures",
    )
    parser.add_argument(
        "--format",
        type=str,
        default="both",
        choices=["png", "svg", "both"],
        help="Output format",
    )
    parser.add_argument(
        "--chr4-alpha",
        type=float,
        default=0.25,
        help="Alpha transparency for chr4 hotspot region (default: 0.25)",
    )
    parser.add_argument(
        "--dpi", type=int, default=300, help="DPI for PNG output (default: 300)"
    )

    return parser.parse_args()


def load_data(input_file, mode):
    """Load QTL overlap summary data and select columns based on mode."""
    print(f"Loading data from {input_file}...")
    df = pd.read_csv(input_file)

    # Select columns based on mode
    if mode == "exact":
        qtl_cols = ["exact_piQTL", "exact_pQTL", "exact_eQTL"]
    else:  # colocal
        qtl_cols = ["colocal_piQTL", "colocal_pQTL", "colocal_eQTL"]

    print(f"  Loaded {len(df)} SNPs")
    print(f"  Using columns: {qtl_cols}")

    # Calculate statistics
    for col in qtl_cols:
        n_nonzero = (df[col] > 0).sum()
        pct_nonzero = 100 * n_nonzero / len(df)
        print(f"  {col}: {n_nonzero} SNPs ({pct_nonzero:.1f}%) with QTLs")

    return df, qtl_cols


def prepare_genomic_positions(df, qtl_cols):
    """
    Create a scaled genomic position coordinate for plotting.
    Chromosome width is proportional to the number of valid SNPs (those with QTLs).
    """
    df = df.copy()

    # Identify valid SNPs (those with at least one QTL > 0)
    df["has_qtl"] = (df[qtl_cols] > 0).any(axis=1)

    # Get chromosome info
    chromosomes = sorted(df["chromosome"].unique())

    # Count valid SNPs per chromosome
    valid_snp_counts = {}
    for chrom in chromosomes:
        chr_mask = df["chromosome"] == chrom
        valid_snp_counts[chrom] = df[chr_mask & df["has_qtl"]].shape[0]

    total_valid_snps = sum(valid_snp_counts.values())
    print(f"  Total valid SNPs (with QTLs): {total_valid_snps}")

    # Create cumulative position (scaled by valid SNP count per chromosome)
    df["genome_pos"] = 0.0
    chr_boundaries = {}
    chr_centers = {}

    cumulative = 0
    for i, chrom in enumerate(chromosomes):
        chr_mask = df["chromosome"] == chrom
        chr_data = df[chr_mask]

        # Width proportional to valid SNP count
        # Use maximum of 1 to avoid division by zero for chromosomes with no valid SNPs
        chr_width = max(1, valid_snp_counts[chrom])

        # Normalize positions within chromosome based on width
        pos_min = chr_data["position"].min()
        pos_max = chr_data["position"].max()
        pos_range = pos_max - pos_min if pos_max > pos_min else 1

        normalized = (chr_data["position"] - pos_min) / pos_range * chr_width
        df.loc[chr_mask, "genome_pos"] = cumulative + normalized

        chr_boundaries[chrom] = (cumulative, cumulative + chr_width)
        chr_centers[chrom] = cumulative + chr_width / 2
        cumulative += chr_width

    return df, chr_boundaries, chr_centers, total_valid_snps


def create_heatmap_data(df, qtl_cols, chr_boundaries, n_bins=200):
    """
    Create heatmap data by binning genome into windows.
    Returns 2D array: 3 QTL types × bins across entire genome.
    """
    total_width = max(bound[1] for bound in chr_boundaries.values())

    # Create matrix: rows=3 QTL types, cols=bins across genome
    heatmap_data = np.zeros((3, n_bins))

    # Bin genome positions
    genome_positions = df["genome_pos"].values
    bin_edges = np.linspace(0, total_width, n_bins + 1)

    # For each QTL type (row)
    for i, qtl_col in enumerate(qtl_cols):
        # Aggregate QTL counts per bin
        for j in range(n_bins):
            bin_mask = (genome_positions >= bin_edges[j]) & (
                genome_positions < bin_edges[j + 1]
            )
            if bin_mask.sum() > 0:
                # Sum QTL counts in this bin for this QTL type
                bin_sum = df.loc[df.index[bin_mask], qtl_col].sum()
                heatmap_data[i, j] = bin_sum

    return heatmap_data


def plot_heatmap(ax, heatmap_data, chr_boundaries, chr_centers, mode, qtl_cols):
    """Plot genome-wide heatmap of QTL counts with per-row normalization."""
    # Apply log transformation for visualization
    heatmap_log = np.log10(heatmap_data + 1)

    # Normalize each row independently for better visualization
    heatmap_normalized = np.zeros_like(heatmap_log)
    for i in range(heatmap_log.shape[0]):
        row = heatmap_log[i, :]
        row_min = row.min()
        row_max = row.max()
        if row_max > row_min:
            heatmap_normalized[i, :] = (row - row_min) / (row_max - row_min)
        else:
            heatmap_normalized[i, :] = row

    # Create heatmap
    chromosomes = sorted(chr_boundaries.keys())
    total_width = max(bound[1] for bound in chr_boundaries.values())

    im = ax.imshow(
        heatmap_normalized,
        aspect="auto",
        cmap="viridis",
        interpolation="nearest",
        extent=[0, total_width, 3, 0],
        vmin=0,
        vmax=1,
    )

    # Set y-axis (QTL types)
    qtl_labels = [col.split("_")[1] for col in qtl_cols]  # Extract piQTL, pQTL, eQTL
    ax.set_yticks([0.5, 1.5, 2.5])
    ax.set_yticklabels(qtl_labels, fontsize=6)
    ax.set_ylabel("QTL type", fontsize=6)

    # Set x-axis to match scatter plots
    ax.set_xlim(0, total_width)
    ax.set_xticks([chr_centers[c] for c in chromosomes])
    ax.set_xticklabels([str(c) for c in chromosomes], fontsize=5)
    ax.set_xlabel("Chromosome", fontsize=6)

    # Add chromosome boundaries
    for chrom in chromosomes[:-1]:
        boundary = chr_boundaries[chrom][1]
        ax.axvline(boundary, color="white", linewidth=0.5, alpha=0.7)

    # Title
    ax.set_title(
        f"QTL Overlap Heatmap ({mode.capitalize()} mode)", fontsize=7, weight="bold"
    )

    # Colorbar in upper right corner
    cax = ax.inset_axes(
        [0.75, 1.5, 0.2, 0.03]
    )  # [x, y, width, height] in axes coordinates
    cbar = plt.colorbar(im, cax=cax, orientation="horizontal")
    cbar.set_label("Normalized intensity (per row)", fontsize=5)
    cbar.ax.tick_params(labelsize=4)

    # Add panel label
    ax.text(-0.08, 1.05, "A", transform=ax.transAxes, fontsize=8, weight="bold")

    return ax


def plot_albert_indicator(ax, df, chr_boundaries, chr_centers):
    """Plot horizontal strip indicating Albert SNPs."""
    chromosomes = sorted(chr_boundaries.keys())
    total_width = max(bound[1] for bound in chr_boundaries.values())

    # Create color array for each SNP
    colors = df["isAlbert_SNP"].map({True: "gold", False: "lightgray"})

    # Plot as scatter with very small height
    for chrom in chromosomes:
        chr_data = df[df["chromosome"] == chrom]
        chr_colors = colors[chr_data.index]

        # Plot within chromosome boundaries
        x_pos = chr_data["genome_pos"].values
        y_pos = np.ones(len(chr_data)) * 0.5

        ax.scatter(
            x_pos, y_pos, c=chr_colors, s=10, marker="|", alpha=0.8, linewidths=0.5
        )

    # Styling
    ax.set_xlim(0, total_width)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_ylabel("Albert", fontsize=5, rotation=0, ha="right", va="center")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    # Add chromosome boundaries
    for chrom in chromosomes[:-1]:
        boundary = chr_boundaries[chrom][1]
        ax.axvline(boundary, color="white", linewidth=0.5, alpha=0.5)

    # Add panel label
    ax.text(
        -0.08, 0.5, "B", transform=ax.transAxes, fontsize=8, weight="bold", va="center"
    )

    return ax


def plot_scatter(ax, df, qtl_col, chr_boundaries, chr_centers, chr4_alpha, panel_label):
    """Plot scatter plot of QTL counts vs genome position with log y-axis."""
    # Filter to SNPs with at least 1 QTL
    plot_data = df[df[qtl_col] >= 1].copy()

    chromosomes = sorted(chr_boundaries.keys())
    total_width = max(bound[1] for bound in chr_boundaries.values())

    if len(plot_data) == 0:
        ax.text(
            0.5,
            0.5,
            "No QTLs detected",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=6,
            style="italic",
        )
        ax.set_xlim(0, total_width)
        ax.set_yscale("log")
        ax.set_ylim(1, 10)
        return ax

    print(f"  Plotting {len(plot_data)} SNPs with {qtl_col} >= 1")

    # Prepare colors based on count values
    counts = plot_data[qtl_col].values

    # Plot with different alpha for chr4
    for chrom in chromosomes:
        chr_data = plot_data[plot_data["chromosome"] == chrom]
        if len(chr_data) == 0:
            continue

        # Adjust alpha for chr4 hotspot
        alpha = chr4_alpha if chrom == 4 else 0.4

        ax.scatter(
            chr_data["genome_pos"],
            chr_data[qtl_col],
            c=chr_data[qtl_col],
            s=1,
            alpha=alpha,
            cmap="viridis",
            vmin=counts.min(),
            vmax=counts.max(),
            rasterized=True,
        )

    # Set log scale for y-axis
    ax.set_yscale("log")

    # Determine y-axis limits
    y_min = max(0.5, counts.min() * 0.8)
    y_max = counts.max() * 1.5
    ax.set_ylim(y_min, y_max)

    # Format y-axis labels as decimal notation (1, 10, 100) instead of exponential (10^0, 10^1, 10^2)
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.ticklabel_format(style="plain", axis="y")

    # X-axis
    ax.set_xlim(0, total_width)
    ax.set_xticks([chr_centers[c] for c in chromosomes])
    ax.set_xticklabels([str(c) for c in chromosomes], fontsize=5)

    # Add chromosome boundaries
    for chrom in chromosomes[:-1]:
        boundary = chr_boundaries[chrom][1]
        ax.axvline(boundary, color="gray", linewidth=0.5, alpha=0.3, linestyle="--")

    # Labels
    qtl_type = qtl_col.split("_")[1]  # Extract piQTL, pQTL, or eQTL
    ax.set_ylabel(f"{qtl_type} count", fontsize=6)
    ax.tick_params(axis="y", labelsize=5)

    # Grid
    ax.grid(axis="y", alpha=0.2, linewidth=0.5)

    # Add panel label
    ax.text(-0.08, 1.05, panel_label, transform=ax.transAxes, fontsize=8, weight="bold")

    return ax


def create_figure(df, qtl_cols, chr_boundaries, chr_centers, mode, chr4_alpha):
    """Create the complete multi-panel figure."""
    print("\nCreating figure...")

    # Figure setup
    fig = plt.figure(figsize=(18 * CM, 16 * CM))

    # Create grid with custom height ratios
    gs = gridspec.GridSpec(
        5, 1, figure=fig, height_ratios=[1, 0.3, 1.5, 1.5, 1.5], hspace=0.4
    )

    # Panel 1: Heatmap
    ax_heatmap = fig.add_subplot(gs[0])
    heatmap_data = create_heatmap_data(df, qtl_cols, chr_boundaries, n_bins=600)
    plot_heatmap(ax_heatmap, heatmap_data, chr_boundaries, chr_centers, mode, qtl_cols)

    # Panel 2: Albert indicator
    ax_albert = fig.add_subplot(gs[1])
    plot_albert_indicator(ax_albert, df, chr_boundaries, chr_centers)

    # Panels 3-5: Scatter plots
    scatter_axes = [fig.add_subplot(gs[i]) for i in range(2, 5)]
    panel_labels = ["C", "D", "E"]

    for i, (ax, qtl_col, label) in enumerate(zip(scatter_axes, qtl_cols, panel_labels)):
        plot_scatter(ax, df, qtl_col, chr_boundaries, chr_centers, chr4_alpha, label)

        # Only show x-axis label on bottom plot
        if i == 2:
            ax.set_xlabel("Chromosome", fontsize=6)
        else:
            ax.set_xticklabels([])

    return fig


def main():
    """Main execution function."""
    args = parse_args()

    print("=" * 80)
    print("QTL Overlap Visualization")
    print("=" * 80)
    print(f"Mode: {args.mode}")
    print(f"Input: {args.input}")
    print(f"Output directory: {args.output_dir}")
    print()

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    df, qtl_cols = load_data(args.input, args.mode)

    # Prepare genomic coordinates
    print("\nPreparing genomic coordinates...")
    df, chr_boundaries, chr_centers, total_valid_snps = prepare_genomic_positions(
        df, qtl_cols
    )
    print(f"  Spanning {len(chr_boundaries)} chromosomes")
    print(
        f"  Total genome width: {max(bound[1] for bound in chr_boundaries.values()):.0f} SNPs"
    )

    # Create figure
    fig = create_figure(
        df, qtl_cols, chr_boundaries, chr_centers, args.mode, args.chr4_alpha
    )

    # Save outputs
    print("\nSaving figures...")
    base_name = f"QTL_overlap_{args.mode}"

    if args.format in ["png", "both"]:
        png_path = output_dir / f"{base_name}.png"
        fig.savefig(png_path, dpi=args.dpi, bbox_inches="tight")
        print(f"  Saved PNG: {png_path}")

    if args.format in ["svg", "both"]:
        svg_path = output_dir / f"{base_name}.svg"
        fig.savefig(svg_path, format="svg", bbox_inches="tight")
        print(f"  Saved SVG: {svg_path}")

    plt.close(fig)

    print("\n" + "=" * 80)
    print("Completed successfully!")
    print("=" * 80)


if __name__ == "__main__":
    main()
