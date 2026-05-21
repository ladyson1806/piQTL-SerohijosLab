"""
Calculate empirical p-values for QTL overlap enrichment.

Compares actual observed overlap counts against the distribution from randomized simulations.
Empirical p-value = (number of random samples >= observed + 1) / (total simulations + 1)
The +1 adjustment avoids p-value = 0 for finite sampling.
"""

import pandas as pd
import argparse
import os


def calculate_empirical_pvalue(randomized_counts, observed_count):
    """
    Calculate one-tailed empirical p-value for enrichment.
    Uses (n_extreme + 1) / (n_total + 1) to avoid p-value = 0.

    Parameters:
    - randomized_counts: array of counts from randomized simulations
    - observed_count: actual observed count

    Returns:
    - p-value: proportion of random samples >= observed
    """
    n_total = len(randomized_counts)
    n_extreme = sum(randomized_counts >= observed_count)
    p_value = (n_extreme + 1) / (n_total + 1)
    return p_value


def calculate_all_pvalues(overlap_mode, randomized_csv, actual_csv):
    """
    Calculate empirical p-values for all overlap categories.

    Parameters:
    - overlap_mode: 'exact' or 'colocal'
    - randomized_csv: path to randomized overlap counts CSV
    - actual_csv: path to actual overlap status CSV

    Returns:
    - DataFrame with p-values and statistics
    """
    # Read data
    random_df = pd.read_csv(randomized_csv)
    actual_df = pd.read_csv(actual_csv)

    # Extract actual values for this overlap mode
    actual_values = actual_df[actual_df["overlap_mode"] == overlap_mode].iloc[0]

    # Calculate p-values for each category
    results = []

    # pQTL overlap
    pqtl_pvalue = calculate_empirical_pvalue(
        random_df["overlap_pqtl"].values, actual_values["with_pQTL"]
    )
    pqtl_mean = random_df["overlap_pqtl"].mean()
    pqtl_std = random_df["overlap_pqtl"].std()

    results.append(
        {
            "overlap_mode": overlap_mode,
            "category": "pQTL",
            "observed": actual_values["with_pQTL"],
            "random_mean": pqtl_mean,
            "random_std": pqtl_std,
            "fold_enrichment": actual_values["with_pQTL"] / pqtl_mean,
            "empirical_pvalue": pqtl_pvalue,
            "n_simulations": len(random_df),
        }
    )

    # eQTL overlap
    eqtl_pvalue = calculate_empirical_pvalue(
        random_df["overlap_eqtl"].values, actual_values["with_eQTL"]
    )
    eqtl_mean = random_df["overlap_eqtl"].mean()
    eqtl_std = random_df["overlap_eqtl"].std()

    results.append(
        {
            "overlap_mode": overlap_mode,
            "category": "eQTL",
            "observed": actual_values["with_eQTL"],
            "random_mean": eqtl_mean,
            "random_std": eqtl_std,
            "fold_enrichment": actual_values["with_eQTL"] / eqtl_mean,
            "empirical_pvalue": eqtl_pvalue,
            "n_simulations": len(random_df),
        }
    )

    # Both (pQTL + eQTL) overlap
    both_pvalue = calculate_empirical_pvalue(
        random_df["overlap_both"].values, actual_values["both"]
    )
    both_mean = random_df["overlap_both"].mean()
    both_std = random_df["overlap_both"].std()

    results.append(
        {
            "overlap_mode": overlap_mode,
            "category": "both",
            "observed": actual_values["both"],
            "random_mean": both_mean,
            "random_std": both_std,
            "fold_enrichment": actual_values["both"] / both_mean,
            "empirical_pvalue": both_pvalue,
            "n_simulations": len(random_df),
        }
    )

    return pd.DataFrame(results)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate empirical p-values for QTL overlap enrichment.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python calculate_empirical_pvalues.py \\
    --randomized-exact out/randomized_overlap_counts_exact/overlap_counts_results.csv \\
    --randomized-colocal out/randomized_overlap_counts_colocal/overlap_counts_results.csv \\
    --actual data/actual_piQTL_SNPs_overlap_status.csv \\
    --output out/empirical_pvalues.csv
        """,
    )

    parser.add_argument(
        "--randomized-exact",
        required=True,
        help="Path to randomized exact overlap counts CSV",
    )
    parser.add_argument(
        "--randomized-colocal",
        required=True,
        help="Path to randomized colocal overlap counts CSV",
    )
    parser.add_argument(
        "--actual", required=True, help="Path to actual overlap status CSV"
    )
    parser.add_argument(
        "--output", required=True, help="Output CSV file path for p-values"
    )

    args = parser.parse_args()

    # Validate input files
    if not os.path.exists(args.randomized_exact):
        parser.error(f"Randomized exact CSV not found: {args.randomized_exact}")
    if not os.path.exists(args.randomized_colocal):
        parser.error(f"Randomized colocal CSV not found: {args.randomized_colocal}")
    if not os.path.exists(args.actual):
        parser.error(f"Actual CSV not found: {args.actual}")

    # Calculate p-values for both overlap modes
    exact_results = calculate_all_pvalues("exact", args.randomized_exact, args.actual)
    colocal_results = calculate_all_pvalues(
        "colocal", args.randomized_colocal, args.actual
    )

    # Combine results
    all_results = pd.concat([exact_results, colocal_results], ignore_index=True)

    # Save to CSV
    all_results.to_csv(args.output, index=False, float_format="%.6f")

    print(f"Empirical p-values saved to: {args.output}")
    print("\nResults summary:")
    print("=" * 100)

    # Display formatted results
    for _, row in all_results.iterrows():
        print(f"\n{row['overlap_mode'].upper()} overlap - {row['category']}:")
        print(f"  Observed count:        {row['observed']:.0f}")
        print(
            f"  Random mean ± SD:      {row['random_mean']:.2f} ± {row['random_std']:.2f}"
        )
        print(f"  Fold enrichment:       {row['fold_enrichment']:.3f}x")
        print(f"  Empirical p-value:     {row['empirical_pvalue']:.6f}")
        print(f"  Number of simulations: {row['n_simulations']:.0f}")

        # Significance indicator
        if row["empirical_pvalue"] < 0.001:
            sig = "***"
        elif row["empirical_pvalue"] < 0.01:
            sig = "**"
        elif row["empirical_pvalue"] < 0.05:
            sig = "*"
        else:
            sig = "ns"
        print(f"  Significance:          {sig}")

    print("\n" + "=" * 100)
    print(
        "Significance levels: *** p < 0.001, ** p < 0.01, * p < 0.05, ns = not significant"
    )
