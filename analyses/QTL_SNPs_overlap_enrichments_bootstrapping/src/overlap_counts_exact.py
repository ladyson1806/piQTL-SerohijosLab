#!/usr/bin/env python3
"""
Count overlapping piQTL SNPs with pQTL and eQTL using random sampling.

This script randomly samples piQTL SNPs multiple times and counts how many
overlap with pQTL, eQTL, and both datasets. Results are saved to CSV files
with per-iteration counts and the randomly selected SNP IDs.
"""

import pandas as pd
import numpy as np
import argparse
from pathlib import Path
from datetime import datetime


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Count overlapping piQTL SNPs with pQTL and eQTL using random sampling."
    )
    parser.add_argument(
        "--n",
        type=int,
        required=True,
        help="Number of random sampling iterations",
    )
    parser.add_argument(
        "--seed",
        type=int,
        required=True,
        help="Random seed for reproducibility",
    )
    parser.add_argument(
        "--num_snp",
        type=int,
        required=True,
        help="Number of piQTL SNPs to sample per iteration",
    )
    parser.add_argument(
        "--piqtl",
        type=str,
        required=True,
        help="Path to piQTL SNP annotation CSV file",
    )
    parser.add_argument(
        "--pqtl",
        type=str,
        required=True,
        help="Path to pQTL results CSV file",
    )
    parser.add_argument(
        "--eqtl",
        type=str,
        required=True,
        help="Path to eQTL results CSV file",
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default="out",
        help="Output directory (default: out)",
    )
    return parser.parse_args()


def load_qtl_data(piqtl_path, pqtl_path, eqtl_path):
    """
    Load QTL data from CSV files.

    Args:
        piqtl_path: Path to piQTL SNP annotation file
        pqtl_path: Path to pQTL results file
        eqtl_path: Path to eQTL results file

    Returns:
        Tuple of (piqtl_df, pqtl_markers, eqtl_markers)
    """
    print("Loading QTL data files...")

    # Load piQTL data (has both SNP ID and SNP_marker)
    piqtl_df = pd.read_csv(piqtl_path)
    print(f"  - piQTL: {len(piqtl_df)} SNPs")

    # Load pQTL and eQTL data (only SNP_marker needed)
    pqtl_df = pd.read_csv(pqtl_path)
    eqtl_df = pd.read_csv(eqtl_path)

    # Extract sets of SNP markers for fast lookup
    pqtl_markers = set(pqtl_df["SNP_marker"].unique())
    eqtl_markers = set(eqtl_df["SNP_marker"].unique())

    print(f"  - pQTL: {len(pqtl_markers)} unique SNP markers")
    print(f"  - eQTL: {len(eqtl_markers)} unique SNP markers")

    return piqtl_df, pqtl_markers, eqtl_markers


def run_overlap_analysis(
    piqtl_df, pqtl_markers, eqtl_markers, n_iterations, num_snp, seed
):
    """
    Run random sampling iterations and count overlaps.

    Args:
        piqtl_df: DataFrame with piQTL data (columns: SNP, SNP_marker, ...)
        pqtl_markers: Set of pQTL SNP markers
        eqtl_markers: Set of eQTL SNP markers
        n_iterations: Number of random sampling iterations
        num_snp: Number of piQTL SNPs to sample per iteration
        seed: Random seed for reproducibility

    Returns:
        Tuple of (results_list, selected_snps_list)
    """
    np.random.seed(seed)

    # Create mapping from SNP_marker to SNP ID
    marker_to_id = dict(zip(piqtl_df["SNP_marker"], piqtl_df["SNP"]))

    # All available piQTL SNP markers
    all_piqtl_markers = piqtl_df["SNP_marker"].values

    if num_snp > len(all_piqtl_markers):
        raise ValueError(
            f"num_snp ({num_snp}) cannot exceed total piQTL SNPs ({len(all_piqtl_markers)})"
        )

    results_list = []
    selected_snps_list = []

    print(f"\nRunning {n_iterations} random sampling iterations...")
    print(f"  - Sampling {num_snp} piQTL SNPs per iteration")
    print(f"  - Random seed: {seed}")

    for iteration in range(1, n_iterations + 1):
        # Randomly sample piQTL SNP markers
        sampled_markers = np.random.choice(
            all_piqtl_markers, size=num_snp, replace=False
        )

        # Convert markers to SNP IDs for output
        sampled_ids = sorted([int(marker_to_id[marker]) for marker in sampled_markers])
        sampled_ids = map(str, sampled_ids)  # Convert to strings for output

        # Count overlaps
        overlap_pqtl = sum(1 for marker in sampled_markers if marker in pqtl_markers)
        overlap_eqtl = sum(1 for marker in sampled_markers if marker in eqtl_markers)
        overlap_both = sum(
            1
            for marker in sampled_markers
            if marker in pqtl_markers and marker in eqtl_markers
        )

        # Store results
        results_list.append(
            {
                "iteration": iteration,
                "overlap_pqtl": overlap_pqtl,
                "overlap_eqtl": overlap_eqtl,
                "overlap_both": overlap_both,
            }
        )
        selected_snps_list.append(",".join(sampled_ids))

        if iteration % max(1, n_iterations // 10) == 0 or iteration == 1:
            print(
                f"  - Iteration {iteration}/{n_iterations}: pQTL={overlap_pqtl}, eQTL={overlap_eqtl}, Both={overlap_both}"
            )

    return results_list, selected_snps_list


def write_results(results_list, selected_snps_list, outdir):
    """
    Write results to CSV files.

    Args:
        results_list: List of dictionaries with overlap counts
        selected_snps_list: List of comma-separated SNP IDs per iteration
        outdir: Output directory
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Write overlap counts
    results_df = pd.DataFrame(results_list)
    overlap_counts_file = outdir / "overlap_counts_results.csv"
    results_df.to_csv(overlap_counts_file, index=False)
    print(f"\n✓ Overlap counts saved to: {overlap_counts_file}")

    # Write randomly selected SNP IDs
    selected_snps_file = outdir / "randomly_selected_piQTL_SNPs.csv"
    with open(selected_snps_file, "w") as f:
        for snps in selected_snps_list:
            f.write(snps + "\n")
    print(f"✓ Randomly selected piQTL SNPs saved to: {selected_snps_file}")

    return overlap_counts_file, selected_snps_file


def write_config_log(
    config_log_file,
    seed,
    n_iterations,
    num_snp,
    piqtl_path,
    pqtl_path,
    eqtl_path,
    piqtl_count,
    pqtl_count,
    eqtl_count,
    results_list,
    overlap_counts_file,
    selected_snps_file,
):
    """
    Write configuration and summary statistics to log file.

    Args:
        config_log_file: Path to output config log file
        seed: Random seed used
        n_iterations: Number of iterations
        num_snp: Number of SNPs sampled per iteration
        piqtl_path, pqtl_path, eqtl_path: Input file paths
        piqtl_count, pqtl_count, eqtl_count: Data counts
        results_list: List of overlap count dictionaries
        overlap_counts_file: Path to overlap counts output
        selected_snps_file: Path to selected SNPs output
    """
    results_df = pd.DataFrame(results_list)

    # Calculate statistics
    stats_pqtl = {
        "mean": results_df["overlap_pqtl"].mean(),
        "std": results_df["overlap_pqtl"].std(),
        "min": results_df["overlap_pqtl"].min(),
        "max": results_df["overlap_pqtl"].max(),
    }
    stats_eqtl = {
        "mean": results_df["overlap_eqtl"].mean(),
        "std": results_df["overlap_eqtl"].std(),
        "min": results_df["overlap_eqtl"].min(),
        "max": results_df["overlap_eqtl"].max(),
    }
    stats_both = {
        "mean": results_df["overlap_both"].mean(),
        "std": results_df["overlap_both"].std(),
        "min": results_df["overlap_both"].min(),
        "max": results_df["overlap_both"].max(),
    }

    with open(config_log_file, "w") as f:
        f.write("=" * 80 + "\n")
        f.write("OVERLAP ANALYSIS CONFIGURATION AND SUMMARY\n")
        f.write("=" * 80 + "\n")
        f.write(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("\n")

        # Configuration
        f.write("CONFIGURATION:\n")
        f.write("-" * 80 + "\n")
        f.write(f"Random seed: {seed}\n")
        f.write(f"Number of iterations: {n_iterations}\n")
        f.write(f"Number of piQTL SNPs sampled per iteration: {num_snp}\n")
        f.write("\n")

        # Input files
        f.write("INPUT FILES:\n")
        f.write("-" * 80 + "\n")
        f.write(f"piQTL annotation file: {piqtl_path}\n")
        f.write(f"  - Total SNPs: {piqtl_count}\n")
        f.write(f"pQTL results file: {pqtl_path}\n")
        f.write(f"  - Total unique SNP markers: {pqtl_count}\n")
        f.write(f"eQTL results file: {eqtl_path}\n")
        f.write(f"  - Total unique SNP markers: {eqtl_count}\n")
        f.write("\n")

        # Output files
        f.write("OUTPUT FILES:\n")
        f.write("-" * 80 + "\n")
        f.write(f"Overlap counts: {overlap_counts_file}\n")
        f.write(f"Selected SNPs: {selected_snps_file}\n")
        f.write(f"Config log: {config_log_file}\n")
        f.write("\n")

        # Summary statistics
        f.write("SUMMARY STATISTICS:\n")
        f.write("-" * 80 + "\n")

        f.write("\npiQTL → pQTL overlaps:\n")
        f.write(f"  Mean: {stats_pqtl['mean']:.2f} ± {stats_pqtl['std']:.2f}\n")
        f.write(f"  Range: {stats_pqtl['min']} - {stats_pqtl['max']}\n")
        f.write(
            f"  Percentage of sampled SNPs: {100 * stats_pqtl['mean'] / num_snp:.2f}% ± {100 * stats_pqtl['std'] / num_snp:.2f}%\n"
        )

        f.write("\npiQTL → eQTL overlaps:\n")
        f.write(f"  Mean: {stats_eqtl['mean']:.2f} ± {stats_eqtl['std']:.2f}\n")
        f.write(f"  Range: {stats_eqtl['min']} - {stats_eqtl['max']}\n")
        f.write(
            f"  Percentage of sampled SNPs: {100 * stats_eqtl['mean'] / num_snp:.2f}% ± {100 * stats_eqtl['std'] / num_snp:.2f}%\n"
        )

        f.write("\npiQTL overlaps with BOTH pQTL AND eQTL:\n")
        f.write(f"  Mean: {stats_both['mean']:.2f} ± {stats_both['std']:.2f}\n")
        f.write(f"  Range: {stats_both['min']} - {stats_both['max']}\n")
        f.write(
            f"  Percentage of sampled SNPs: {100 * stats_both['mean'] / num_snp:.2f}% ± {100 * stats_both['std'] / num_snp:.2f}%\n"
        )

        f.write("\n" + "=" * 80 + "\n")


def main():
    """Main function."""
    args = parse_arguments()

    # Load data
    piqtl_df, pqtl_markers, eqtl_markers = load_qtl_data(
        args.piqtl, args.pqtl, args.eqtl
    )

    # Run overlap analysis
    results_list, selected_snps_list = run_overlap_analysis(
        piqtl_df, pqtl_markers, eqtl_markers, args.n, args.num_snp, args.seed
    )

    # Write results
    overlap_counts_file, selected_snps_file = write_results(
        results_list, selected_snps_list, args.outdir
    )

    # Write config log
    config_log_file = Path(args.outdir) / "config.txt"
    write_config_log(
        config_log_file,
        args.seed,
        args.n,
        args.num_snp,
        args.piqtl,
        args.pqtl,
        args.eqtl,
        len(piqtl_df),
        len(pqtl_markers),
        len(eqtl_markers),
        results_list,
        overlap_counts_file,
        selected_snps_file,
    )
    print(f"✓ Config log saved to: {config_log_file}\n")

    print("✓ Analysis complete!")


if __name__ == "__main__":
    main()
