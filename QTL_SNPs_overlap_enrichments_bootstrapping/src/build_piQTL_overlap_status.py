#!/usr/bin/env python3
"""
Build a per-SNP overlap status table for all piQTL SNPs.
"""

import argparse
from pathlib import Path

import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Build overlap status table for all piQTL SNPs."
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
        "--actual",
        type=str,
        required=True,
        help="Path to actual piQTL SNP ID list (one per line)",
    )
    parser.add_argument(
        "--out",
        type=str,
        required=True,
        help="Output CSV path",
    )
    return parser.parse_args()


def is_colocal(pos_a, start_a, end_a, pos_b, start_b, end_b):
    return start_b <= pos_a <= end_b and start_a <= pos_b <= end_a


def build_qtl_by_chr(qtl_df):
    qtl_by_chr = {}
    for _, row in qtl_df.iterrows():
        chr_val = int(row["chr"])
        qtl_by_chr.setdefault(chr_val, []).append(
            (int(row["pos"]), int(row["start"]), int(row["end"]))
        )
    return qtl_by_chr


def load_actual_ids(actual_path):
    actual_ids = set()
    with open(actual_path, "r") as handle:
        for line in handle:
            value = line.strip()
            if value:
                actual_ids.add(int(value))
    return actual_ids


def build_overlap_status(piqtl_df, pqtl_df, eqtl_df, actual_ids):
    pqtl_markers = set(pqtl_df["SNP_marker"].unique())
    eqtl_markers = set(eqtl_df["SNP_marker"].unique())

    pqtl_by_chr = build_qtl_by_chr(pqtl_df)
    eqtl_by_chr = build_qtl_by_chr(eqtl_df)

    rows = []

    for _, row in piqtl_df.iterrows():
        snp_id = int(row["SNP"])
        snp_marker = row["SNP_marker"]
        chr_val = int(row["chr"])
        pos = int(row["pos"])
        start = int(row["start"])
        end = int(row["end"])

        exact_pqtl = snp_marker in pqtl_markers
        exact_eqtl = snp_marker in eqtl_markers
        exact_both = exact_pqtl and exact_eqtl

        pqtl_rows = pqtl_by_chr.get(chr_val, [])
        eqtl_rows = eqtl_by_chr.get(chr_val, [])

        colocal_pqtl = any(
            is_colocal(pos, start, end, p_pos, p_start, p_end)
            for p_pos, p_start, p_end in pqtl_rows
        )
        colocal_eqtl = any(
            is_colocal(pos, start, end, e_pos, e_start, e_end)
            for e_pos, e_start, e_end in eqtl_rows
        )
        colocal_both = colocal_pqtl and colocal_eqtl

        rows.append(
            {
                "SNP": snp_id,
                "SNP_marker": snp_marker,
                "exact_pQTL_overlap": exact_pqtl,
                "exact_eQTL_overlap": exact_eqtl,
                "exact_both_overlap": exact_both,
                "colocal_pQTL_overlap": colocal_pqtl,
                "colocal_eQTL_overlap": colocal_eqtl,
                "colocal_both_overlap": colocal_both,
                "Is_actual_piQTL": snp_id in actual_ids,
            }
        )

    return pd.DataFrame(rows)


def main():
    args = parse_arguments()

    piqtl_df = pd.read_csv(args.piqtl)
    pqtl_df = pd.read_csv(args.pqtl)
    eqtl_df = pd.read_csv(args.eqtl)
    actual_ids = load_actual_ids(args.actual)

    result_df = build_overlap_status(piqtl_df, pqtl_df, eqtl_df, actual_ids)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    result_df.to_csv(out_path, index=False)
    print(f"Wrote overlap status table: {out_path}")


if __name__ == "__main__":
    main()
