#!/usr/bin/env python3
"""
Count overlapping piQTL SNPs with pQTL and eQTL using co-local overlap.

This script randomly samples piQTL SNPs multiple times and counts how many
are co-localized with pQTL, eQTL, and both datasets. Co-local overlap requires
mutual containment of peak positions within each other's ranges.
"""

import argparse
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Count overlapping piQTL SNPs with pQTL and eQTL using co-local overlap."
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


def is_colocal(pos_a, start_a, end_a, pos_b, start_b, end_b):
    """Return True if peaks are mutually contained within each other's ranges."""
    return start_b <= pos_a <= end_b and start_a <= pos_b <= end_a


def build_qtl_by_chr(qtl_df):
    """Group QTL rows by chromosome for faster lookup."""
    qtl_by_chr = {}
    for _, row in qtl_df.iterrows():
        chr_val = int(row["chr"])
        qtl_by_chr.setdefault(chr_val, []).append(
            (int(row["pos"]), int(row["start"]), int(row["end"]))
        )
    return qtl_by_chr


def load_qtl_data(piqtl_path, pqtl_path, eqtl_path):
    """Load QTL data from CSV files."""
    print("Loading QTL data files...")

    piqtl_df = pd.read_csv(piqtl_path)
    pqtl_df = pd.read_csv(pqtl_path)
    eqtl_df = pd.read_csv(eqtl_path)

    piqtl_df = piqtl_df.drop_duplicates(subset=["SNP_marker"])

    print(f"  - piQTL: {len(piqtl_df)} SNPs")
    print(f"  - pQTL: {len(pqtl_df)} rows")
    print(f"  - eQTL: {len(eqtl_df)} rows")

    pqtl_by_chr = build_qtl_by_chr(pqtl_df)
    eqtl_by_chr = build_qtl_by_chr(eqtl_df)

    return piqtl_df, pqtl_by_chr, eqtl_by_chr


def run_overlap_analysis(
    piqtl_df, pqtl_by_chr, eqtl_by_chr, n_iterations, num_snp, seed
):
    """Run random sampling iterations and count co-local overlaps."""
    np.random.seed(seed)

    piqtl_map = {}
    for _, row in piqtl_df.iterrows():
        piqtl_map[row["SNP_marker"]] = (
            int(row["SNP"]),
            int(row["chr"]),
            int(row["pos"]),
            int(row["start"]),
            int(row["end"]),
        )

    all_piqtl_markers = np.array(list(piqtl_map.keys()))

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
        sampled_markers = np.random.choice(
            all_piqtl_markers, size=num_snp, replace=False
        )

        overlap_pqtl = 0
        overlap_eqtl = 0
        overlap_both = 0

        sampled_ids = []

        for marker in sampled_markers:
            snp_id, chr_val, pos, start, end = piqtl_map[marker]
            sampled_ids.append(snp_id)

            pqtl_rows = pqtl_by_chr.get(chr_val, [])
            eqtl_rows = eqtl_by_chr.get(chr_val, [])

            is_pqtl = any(
                is_colocal(pos, start, end, p_pos, p_start, p_end)
                for p_pos, p_start, p_end in pqtl_rows
            )
            is_eqtl = any(
                is_colocal(pos, start, end, e_pos, e_start, e_end)
                for e_pos, e_start, e_end in eqtl_rows
            )

            if is_pqtl:
                overlap_pqtl += 1
            if is_eqtl:
                overlap_eqtl += 1
            if is_pqtl and is_eqtl:
                overlap_both += 1

        sampled_ids = map(str, sorted(sampled_ids))

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
    """Write results to CSV files."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    results_df = pd.DataFrame(results_list)
    overlap_counts_file = outdir / "overlap_counts_results.csv"
    results_df.to_csv(overlap_counts_file, index=False)
    print(f"\n✓ Overlap counts saved to: {overlap_counts_file}")

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
    pqtl_rows,
    eqtl_rows,
    results_list,
    overlap_counts_file,
    selected_snps_file,
):
    """Write configuration and summary statistics to log file."""
    results_df = pd.DataFrame(results_list)

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
        f.write("CO-LOCAL OVERLAP ANALYSIS CONFIGURATION AND SUMMARY\n")
        f.write("=" * 80 + "\n")
        f.write(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("\n")

        f.write("CONFIGURATION:\n")
        f.write("-" * 80 + "\n")
        f.write(f"Random seed: {seed}\n")
        f.write(f"Number of iterations: {n_iterations}\n")
        f.write(f"Number of piQTL SNPs sampled per iteration: {num_snp}\n")
        f.write("\n")

        f.write("INPUT FILES:\n")
        f.write("-" * 80 + "\n")
        f.write(f"piQTL annotation file: {piqtl_path}\n")
        f.write(f"  - Total SNPs: {piqtl_count}\n")
        f.write(f"pQTL results file: {pqtl_path}\n")
        f.write(f"  - Total rows: {pqtl_rows}\n")
        f.write(f"eQTL results file: {eqtl_path}\n")
        f.write(f"  - Total rows: {eqtl_rows}\n")
        f.write("\n")

        f.write("OUTPUT FILES:\n")
        f.write("-" * 80 + "\n")
        f.write(f"Overlap counts: {overlap_counts_file}\n")
        f.write(f"Selected SNPs: {selected_snps_file}\n")
        f.write(f"Config log: {config_log_file}\n")
        f.write("\n")

        f.write("SUMMARY STATISTICS:\n")
        f.write("-" * 80 + "\n")

        f.write("\npiQTL → pQTL co-local overlaps:\n")
        f.write(f"  Mean: {stats_pqtl['mean']:.2f} ± {stats_pqtl['std']:.2f}\n")
        f.write(f"  Range: {stats_pqtl['min']} - {stats_pqtl['max']}\n")
        f.write(
            f"  Percentage of sampled SNPs: {100 * stats_pqtl['mean'] / num_snp:.2f}% ± {100 * stats_pqtl['std'] / num_snp:.2f}%\n"
        )

        f.write("\npiQTL → eQTL co-local overlaps:\n")
        f.write(f"  Mean: {stats_eqtl['mean']:.2f} ± {stats_eqtl['std']:.2f}\n")
        f.write(f"  Range: {stats_eqtl['min']} - {stats_eqtl['max']}\n")
        f.write(
            f"  Percentage of sampled SNPs: {100 * stats_eqtl['mean'] / num_snp:.2f}% ± {100 * stats_eqtl['std'] / num_snp:.2f}%\n"
        )

        f.write("\npiQTL co-local overlaps with BOTH pQTL AND eQTL:\n")
        f.write(f"  Mean: {stats_both['mean']:.2f} ± {stats_both['std']:.2f}\n")
        f.write(f"  Range: {stats_both['min']} - {stats_both['max']}\n")
        f.write(
            f"  Percentage of sampled SNPs: {100 * stats_both['mean'] / num_snp:.2f}% ± {100 * stats_both['std'] / num_snp:.2f}%\n"
        )

        f.write("\n" + "=" * 80 + "\n")


def main():
    """Main function."""
    args = parse_arguments()

    piqtl_df, pqtl_by_chr, eqtl_by_chr = load_qtl_data(args.piqtl, args.pqtl, args.eqtl)

    results_list, selected_snps_list = run_overlap_analysis(
        piqtl_df, pqtl_by_chr, eqtl_by_chr, args.n, args.num_snp, args.seed
    )

    overlap_counts_file, selected_snps_file = write_results(
        results_list, selected_snps_list, args.outdir
    )

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
        sum(len(v) for v in pqtl_by_chr.values()),
        sum(len(v) for v in eqtl_by_chr.values()),
        results_list,
        overlap_counts_file,
        selected_snps_file,
    )
    print(f"✓ Config log saved to: {config_log_file}\n")

    print("✓ Analysis complete!")


if __name__ == "__main__":
    main()
