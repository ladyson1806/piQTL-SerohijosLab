#!/usr/bin/env python3
"""
Calculate empirical p-values for reference values against bootstrap simulations.

This script outputs two one-tailed empirical p-value tables with +1 correction:

1) GE table (simulation >= reference)
    p_ge = (count(sim >= ref) + 1) / (n + 1)

2) Smaller-than table (reference < simulation)
    p_ref_smaller = (count(sim > ref) + 1) / (n + 1)

For each metric (mean_diff, median_diff, cis_trans_ratio, trans_pct):
- pQTL and eQTL reference values are tested against their respective simulations
- piQTL values are tested against both pQTL and eQTL simulations
- Tests are reported for global and piQTL_target_subset scopes
"""

import argparse
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate empirical p-values for reference values against bootstrap simulations"
    )
    parser.add_argument(
        "--pqtl-sim", required=True, help="pQTL simulation_summary.tsv file"
    )
    parser.add_argument(
        "--eqtl-sim", required=True, help="eQTL simulation_summary.tsv file"
    )
    parser.add_argument(
        "--summary-table", required=True, help="cis_trans_qtl_summary_combined.csv"
    )
    parser.add_argument(
        "--output-pvalues-ge",
        default="out/tables/cis_trans_qtl_pvalues.csv",
        help="Output CSV for p(sim >= ref) with +1 correction",
    )
    parser.add_argument(
        "--output-pvalues-smaller",
        default="out/tables/cis_trans_qtl_pvalues_ref_smaller.csv",
        help="Output CSV for p(ref < sim) with +1 correction",
    )
    return parser.parse_args()


def load_simulation_data(pqtl_path, eqtl_path):
    """Load and prepare simulation data."""
    pqtl_sim = pd.read_csv(pqtl_path, sep="\t")
    eqtl_sim = pd.read_csv(eqtl_path, sep="\t")

    # Calculate trans_pct for simulations
    pqtl_sim["trans_pct"] = (
        100.0
        * pqtl_sim["count_trans"]
        / (pqtl_sim["count_cis"] + pqtl_sim["count_trans"])
    )
    eqtl_sim["trans_pct"] = (
        100.0
        * eqtl_sim["count_trans"]
        / (eqtl_sim["count_cis"] + eqtl_sim["count_trans"])
    )

    return pqtl_sim, eqtl_sim


def load_reference_markers(summary_path):
    """Load reference values for markers."""
    summary_df = pd.read_csv(summary_path)

    markers = {}
    for qtl_type in ["piQTL", "pQTL", "eQTL"]:
        for scope in ["global", "piQTL_target_subset"]:
            row = summary_df[
                (summary_df["qtl_type"] == qtl_type) & (summary_df["scope"] == scope)
            ]
            if not row.empty:
                key = f"{qtl_type}_{scope}"
                markers[key] = {
                    "mean_diff": row["mean_diff_cis_minus_trans"].values[0],
                    "median_diff": row["median_diff_cis_minus_trans"].values[0],
                    "cis_trans_ratio": row["cis_trans_ratio"].values[0],
                    "trans_pct": row["trans_pct"].values[0],
                }

    return markers


def calculate_empirical_pvalues(reference_value, simulations):
    """Calculate two one-tailed empirical p-values with +1 correction."""
    n_total = len(simulations)
    n_greater_equal = int((simulations >= reference_value).sum())
    n_sim_greater = int((simulations <= reference_value).sum())

    p_ge = (n_greater_equal + 1) / (n_total + 1)
    p_ref_smaller = (n_sim_greater + 1) / (n_total + 1)

    return {
        "p_ge": p_ge,
        "p_ref_smaller": p_ref_smaller,
        "n_ge": n_greater_equal,
        "n_sim_gt": n_sim_greater,
        "n_total": n_total,
    }


def reorder_columns(df, prefixes):
    col_order = ["metric", "metric_label", "qtl_type", "scope", "reference_value"]
    for prefix in prefixes:
        for col in df.columns:
            if col not in col_order and col.startswith(prefix):
                col_order.append(col)
    return df[col_order]


def main():
    args = parse_args()

    print("Loading simulation data...")
    pqtl_sim, eqtl_sim = load_simulation_data(args.pqtl_sim, args.eqtl_sim)

    print("Loading reference markers...")
    markers = load_reference_markers(args.summary_table)

    print("\n" + "=" * 80)
    print("Calculating Empirical P-values (+1 correction)")
    print("=" * 80)

    results_ge = []
    results_smaller = []

    metrics = [
        ("mean_diff", "Mean difference (cis - trans)"),
        ("median_diff", "Median difference (cis - trans)"),
        ("cis_trans_ratio", "cis/trans ratio"),
        ("trans_pct", "trans-QTL percentage (%)"),
    ]

    for metric_col, metric_label in metrics:
        print(f"\n{metric_label}:")
        print("-" * 80)

        # Test pQTL values
        for scope in ["global", "piQTL_target_subset"]:
            key = f"pQTL_{scope}"
            if key in markers:
                ref_value = markers[key][metric_col]
                pval_result = calculate_empirical_pvalues(
                    ref_value, pqtl_sim[metric_col]
                )

                results_ge.append(
                    {
                        "metric": metric_col,
                        "metric_label": metric_label,
                        "qtl_type": "pQTL",
                        "scope": scope,
                        "reference_value": ref_value,
                        "p_value": pval_result["p_ge"],
                        "n_ge": pval_result["n_ge"],
                        "n_total": pval_result["n_total"],
                    }
                )

                results_smaller.append(
                    {
                        "metric": metric_col,
                        "metric_label": metric_label,
                        "qtl_type": "pQTL",
                        "scope": scope,
                        "reference_value": ref_value,
                        "p_value_ref_smaller": pval_result["p_ref_smaller"],
                        "n_sim_gt": pval_result["n_sim_gt"],
                        "n_total": pval_result["n_total"],
                    }
                )

                print(
                    f"  pQTL ({scope}):       ref={ref_value:.6f},  p_ge={pval_result['p_ge']:.6f}, p_ref_smaller={pval_result['p_ref_smaller']:.6f}"
                )

        # Test eQTL values
        for scope in ["global", "piQTL_target_subset"]:
            key = f"eQTL_{scope}"
            if key in markers:
                ref_value = markers[key][metric_col]
                pval_result = calculate_empirical_pvalues(
                    ref_value, eqtl_sim[metric_col]
                )

                results_ge.append(
                    {
                        "metric": metric_col,
                        "metric_label": metric_label,
                        "qtl_type": "eQTL",
                        "scope": scope,
                        "reference_value": ref_value,
                        "p_value": pval_result["p_ge"],
                        "n_ge": pval_result["n_ge"],
                        "n_total": pval_result["n_total"],
                    }
                )

                results_smaller.append(
                    {
                        "metric": metric_col,
                        "metric_label": metric_label,
                        "qtl_type": "eQTL",
                        "scope": scope,
                        "reference_value": ref_value,
                        "p_value_ref_smaller": pval_result["p_ref_smaller"],
                        "n_sim_gt": pval_result["n_sim_gt"],
                        "n_total": pval_result["n_total"],
                    }
                )

                print(
                    f"  eQTL ({scope}):       ref={ref_value:.6f},  p_ge={pval_result['p_ge']:.6f}, p_ref_smaller={pval_result['p_ref_smaller']:.6f}"
                )

        # Test piQTL values (against both simulations for reference)
        for scope in ["global", "piQTL_target_subset"]:
            key = f"piQTL_{scope}"
            if key in markers:
                ref_value = markers[key][metric_col]

                # piQTL vs pQTL
                pval_pqtl = calculate_empirical_pvalues(ref_value, pqtl_sim[metric_col])
                pval_eqtl = calculate_empirical_pvalues(ref_value, eqtl_sim[metric_col])

                results_ge.append(
                    {
                        "metric": metric_col,
                        "metric_label": metric_label,
                        "qtl_type": "piQTL",
                        "scope": scope,
                        "reference_value": ref_value,
                        "p_value_vs_pqtl": pval_pqtl["p_ge"],
                        "n_ge_vs_pqtl": pval_pqtl["n_ge"],
                        "n_total_vs_pqtl": pval_pqtl["n_total"],
                        "p_value_vs_eqtl": pval_eqtl["p_ge"],
                        "n_ge_vs_eqtl": pval_eqtl["n_ge"],
                        "n_total_vs_eqtl": pval_eqtl["n_total"],
                    }
                )

                results_smaller.append(
                    {
                        "metric": metric_col,
                        "metric_label": metric_label,
                        "qtl_type": "piQTL",
                        "scope": scope,
                        "reference_value": ref_value,
                        "p_value_ref_smaller_vs_pqtl": pval_pqtl["p_ref_smaller"],
                        "n_sim_gt_vs_pqtl": pval_pqtl["n_sim_gt"],
                        "n_total_vs_pqtl": pval_pqtl["n_total"],
                        "p_value_ref_smaller_vs_eqtl": pval_eqtl["p_ref_smaller"],
                        "n_sim_gt_vs_eqtl": pval_eqtl["n_sim_gt"],
                        "n_total_vs_eqtl": pval_eqtl["n_total"],
                    }
                )

                print(
                    f"  piQTL ({scope}):     ref={ref_value:.6f},  p_ge_vs_pqtl={pval_pqtl['p_ge']:.6f}, p_ge_vs_eqtl={pval_eqtl['p_ge']:.6f}, p_ref_smaller_vs_pqtl={pval_pqtl['p_ref_smaller']:.6f}, p_ref_smaller_vs_eqtl={pval_eqtl['p_ref_smaller']:.6f}"
                )

    ge_df = reorder_columns(pd.DataFrame(results_ge), ["p_value", "n_ge", "n_total"])
    smaller_df = reorder_columns(
        pd.DataFrame(results_smaller),
        ["p_value_ref_smaller", "n_sim_gt", "n_total"],
    )

    output_ge_path = Path(args.output_pvalues_ge)
    output_smaller_path = Path(args.output_pvalues_smaller)
    output_ge_path.parent.mkdir(parents=True, exist_ok=True)
    output_smaller_path.parent.mkdir(parents=True, exist_ok=True)
    ge_df.to_csv(output_ge_path, index=False)
    smaller_df.to_csv(output_smaller_path, index=False)

    print("\n" + "=" * 80)
    print("P-value Calculation Summary (+1 correction)")
    print("=" * 80)
    print(f"\nSaved {len(ge_df)} GE-table rows to: {output_ge_path}")
    print(f"Saved {len(smaller_df)} smaller-than-table rows to: {output_smaller_path}")
    print("\nResults by metric:")

    for metric_col, metric_label in metrics:
        metric_results_ge = ge_df[ge_df["metric"] == metric_col]
        metric_results_smaller = smaller_df[smaller_df["metric"] == metric_col]
        print(f"\n  {metric_label}:")
        for _, row in metric_results_ge.iterrows():
            if row["qtl_type"] == "piQTL":
                print(
                    f"    {row['qtl_type']} ({row['scope']}) [GE]:     p_vs_pqtl={row.get('p_value_vs_pqtl', 'N/A')}, p_vs_eqtl={row.get('p_value_vs_eqtl', 'N/A')}"
                )
            else:
                p_val = row.get("p_value", "N/A")
                print(f"    {row['qtl_type']} ({row['scope']}) [GE]:     p={p_val}")

        for _, row in metric_results_smaller.iterrows():
            if row["qtl_type"] == "piQTL":
                print(
                    f"    {row['qtl_type']} ({row['scope']}) [ref<sim]: p_vs_pqtl={row.get('p_value_ref_smaller_vs_pqtl', 'N/A')}, p_vs_eqtl={row.get('p_value_ref_smaller_vs_eqtl', 'N/A')}"
                )
            else:
                p_val = row.get("p_value_ref_smaller", "N/A")
                print(f"    {row['qtl_type']} ({row['scope']}) [ref<sim]: p={p_val}")

    print("\n" + "=" * 80)
    print("Interpretation of p-values (+1 correction):")
    print("  GE table:          p_ge = (count(sim >= ref) + 1) / (n + 1)")
    print("  Smaller-than table: p_ref_smaller = (count(sim > ref) + 1) / (n + 1)")
    print("=" * 80)


if __name__ == "__main__":
    main()
