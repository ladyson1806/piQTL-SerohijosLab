#!/usr/bin/env python3
"""
Summarize cis/trans effect-size patterns across piQTL, pQTL, and eQTL datasets.

This script creates summary tables for two scopes:
1) Global: all rows in each dataset
2) piQTL-target subset: pQTL/eQTL rows restricted to piQTL target ORFs,
   and piQTL rows restricted to PPIs containing piQTL target gene symbols.

Metrics per dataset and scope:
- mean_diff (mean(cis_effect) - mean(trans_effect))
- median_diff (median(cis_effect) - median(trans_effect))
- cis_trans_ratio (n_cis / n_trans)
- trans_pct (100 * n_trans / (n_cis + n_trans))
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize cis/trans QTL effect metrics for piQTL, pQTL, and eQTL"
    )
    parser.add_argument("--piqtl-file", required=True, help="Annotated piQTL CSV file")
    parser.add_argument("--pqtl-file", required=True, help="Jakobson pQTL TSV file")
    parser.add_argument("--eqtl-file", required=True, help="Albert eQTL TSV file")
    parser.add_argument(
        "--target-file", required=True, help="piQTL target gene CSV file"
    )
    parser.add_argument(
        "--out-global", required=True, help="Output CSV path for global summary"
    )
    parser.add_argument(
        "--out-target",
        required=True,
        help="Output CSV path for piQTL-target subset summary",
    )
    parser.add_argument(
        "--out-combined",
        default="",
        help="Optional output CSV path for combined (global + target) summary",
    )
    return parser.parse_args()


def normalize_cis_trans(series):
    return series.astype(str).str.strip().str.lower()


def validate_columns(df, required_columns, label):
    missing = [column for column in required_columns if column not in df.columns]
    if missing:
        missing_str = ", ".join(missing)
        raise ValueError(f"Missing required columns in {label}: {missing_str}")


def build_summary_row(df, qtl_type, scope):
    cis_values = df.loc[df["cis_trans"] == "cis", "effect_size"].dropna()
    trans_values = df.loc[df["cis_trans"] == "trans", "effect_size"].dropna()

    n_cis = int(cis_values.shape[0])
    n_trans = int(trans_values.shape[0])
    n_total = n_cis + n_trans

    mean_cis = cis_values.mean() if n_cis > 0 else np.nan
    mean_trans = trans_values.mean() if n_trans > 0 else np.nan
    median_cis = cis_values.median() if n_cis > 0 else np.nan
    median_trans = trans_values.median() if n_trans > 0 else np.nan

    mean_diff = mean_cis - mean_trans if n_cis > 0 and n_trans > 0 else np.nan
    median_diff = median_cis - median_trans if n_cis > 0 and n_trans > 0 else np.nan
    cis_trans_ratio = (n_cis / n_trans) if n_trans > 0 else np.nan
    trans_pct = (100.0 * n_trans / n_total) if n_total > 0 else np.nan

    return {
        "qtl_type": qtl_type,
        "scope": scope,
        "n_cis": n_cis,
        "n_trans": n_trans,
        "n_total": n_total,
        "mean_effect_cis": mean_cis,
        "mean_effect_trans": mean_trans,
        "median_effect_cis": median_cis,
        "median_effect_trans": median_trans,
        "mean_diff_cis_minus_trans": mean_diff,
        "median_diff_cis_minus_trans": median_diff,
        "cis_trans_ratio": cis_trans_ratio,
        "trans_pct": trans_pct,
    }


def make_pi_target_mask(piqtl_df, target_gene_symbols):
    def has_target_gene(ppi_value):
        genes = str(ppi_value).split("_")
        return any(gene in target_gene_symbols for gene in genes)

    return piqtl_df["PPI"].apply(has_target_gene)


def main():
    args = parse_args()

    piqtl_path = Path(args.piqtl_file)
    pqtl_path = Path(args.pqtl_file)
    eqtl_path = Path(args.eqtl_file)
    target_path = Path(args.target_file)

    out_global_path = Path(args.out_global)
    out_target_path = Path(args.out_target)
    out_combined_path = Path(args.out_combined) if args.out_combined else None

    piqtl_df = pd.read_csv(piqtl_path)
    pqtl_df = pd.read_csv(pqtl_path, sep="\t")
    eqtl_df = pd.read_csv(eqtl_path, sep="\t")
    target_df = pd.read_csv(target_path)

    validate_columns(piqtl_df, ["PPI", "piQTL_EFFECTSIZE", "cis_trans"], "piQTL")
    validate_columns(pqtl_df, ["protein", "cis_trans"], "pQTL")
    validate_columns(eqtl_df, ["gene", "beta", "cis_trans"], "eQTL")
    validate_columns(target_df, ["Gene", "ORF"], "piQTL target genes")

    if "abs_beta" not in pqtl_df.columns:
        if "beta" in pqtl_df.columns:
            pqtl_df["abs_beta"] = pqtl_df["beta"].abs()
        else:
            raise ValueError("Missing required pQTL effect column: abs_beta")

    piqtl_df = piqtl_df.copy()
    pqtl_df = pqtl_df.copy()
    eqtl_df = eqtl_df.copy()

    piqtl_df["effect_size"] = piqtl_df["piQTL_EFFECTSIZE"].abs()
    pqtl_df["effect_size"] = pqtl_df["abs_beta"]
    eqtl_df["effect_size"] = eqtl_df["beta"].abs()

    piqtl_df["cis_trans"] = normalize_cis_trans(piqtl_df["cis_trans"])
    pqtl_df["cis_trans"] = normalize_cis_trans(pqtl_df["cis_trans"])
    eqtl_df["cis_trans"] = normalize_cis_trans(eqtl_df["cis_trans"])

    piqtl_df = piqtl_df[piqtl_df["cis_trans"].isin(["cis", "trans"])].copy()
    pqtl_df = pqtl_df[pqtl_df["cis_trans"].isin(["cis", "trans"])].copy()
    eqtl_df = eqtl_df[eqtl_df["cis_trans"].isin(["cis", "trans"])].copy()

    target_orfs = set(target_df["ORF"].dropna().astype(str))
    target_gene_symbols = set(target_df["Gene"].dropna().astype(str))

    pi_target_mask = make_pi_target_mask(piqtl_df, target_gene_symbols)
    piqtl_target_df = piqtl_df[pi_target_mask].copy()
    pqtl_target_df = pqtl_df[pqtl_df["protein"].astype(str).isin(target_orfs)].copy()
    eqtl_target_df = eqtl_df[eqtl_df["gene"].astype(str).isin(target_orfs)].copy()

    global_rows = [
        build_summary_row(piqtl_df, "piQTL", "global"),
        build_summary_row(pqtl_df, "pQTL", "global"),
        build_summary_row(eqtl_df, "eQTL", "global"),
    ]
    target_rows = [
        build_summary_row(piqtl_target_df, "piQTL", "piQTL_target_subset"),
        build_summary_row(pqtl_target_df, "pQTL", "piQTL_target_subset"),
        build_summary_row(eqtl_target_df, "eQTL", "piQTL_target_subset"),
    ]

    global_df = pd.DataFrame(global_rows)
    target_df = pd.DataFrame(target_rows)
    combined_df = pd.concat([global_df, target_df], ignore_index=True)

    out_global_path.parent.mkdir(parents=True, exist_ok=True)
    out_target_path.parent.mkdir(parents=True, exist_ok=True)
    global_df.to_csv(out_global_path, index=False)
    target_df.to_csv(out_target_path, index=False)

    if out_combined_path is not None:
        out_combined_path.parent.mkdir(parents=True, exist_ok=True)
        combined_df.to_csv(out_combined_path, index=False)

    print("Saved global summary:", out_global_path)
    print("Saved piQTL-target subset summary:", out_target_path)
    if out_combined_path is not None:
        print("Saved combined summary:", out_combined_path)


if __name__ == "__main__":
    main()
