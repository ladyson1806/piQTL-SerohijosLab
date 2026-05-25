#!/usr/bin/env python3
"""
Compare pQTL and eQTL bootstrap simulation results.

Creates 4 figure panels showing distributions from 1,000 random samplings
of 44 genes/proteins from pQTL and eQTL datasets, with reference markers
for global, piQTL-target subset, and piQTL values.

Each panel has 2 boxplots: left=pQTL, right=eQTL
Output: PNG and SVG formats only
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare pQTL and eQTL simulation results with reference markers"
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
        "--output-dir", required=True, help="Output directory for figures"
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


def create_panel_figure(
    pqtl_sim,
    eqtl_sim,
    markers,
    metric_col,
    y_title,
    panel_title,
    output_path,
):
    """Create a single panel with two boxplots (pQTL left, eQTL right)."""
    fig = make_subplots(
        rows=1,
        cols=2,
        subplot_titles=("pQTL Simulations", "eQTL Simulations"),
        horizontal_spacing=0.12,
        shared_yaxes=True,
    )

    # pQTL boxplot
    fig.add_trace(
        go.Box(
            y=pqtl_sim[metric_col],
            name="pQTL",
            marker_color="#1f77b4",
            boxmean=False,
            showlegend=False,
        ),
        row=1,
        col=1,
    )

    # eQTL boxplot
    fig.add_trace(
        go.Box(
            y=eqtl_sim[metric_col],
            name="eQTL",
            marker_color="#ff7f0e",
            boxmean=False,
            showlegend=False,
        ),
        row=1,
        col=2,
    )

    # Add reference markers for pQTL subplot
    if f"pQTL_global" in markers:
        fig.add_hline(
            y=markers["pQTL_global"][metric_col],
            line_dash="dash",
            line_color="red",
            line_width=2,
            annotation_text="Total pQTL",
            annotation_position="right",
            row=1,
            col=1,
        )
    if f"pQTL_piQTL_target_subset" in markers:
        fig.add_hline(
            y=markers["pQTL_piQTL_target_subset"][metric_col],
            line_dash="dot",
            line_color="green",
            line_width=2,
            annotation_text="44 piQTL-restricted pQTL",
            annotation_position="right",
            row=1,
            col=1,
        )
    if f"piQTL_global" in markers:
        fig.add_hline(
            y=markers["piQTL_global"][metric_col],
            line_dash="solid",
            line_color="purple",
            line_width=2,
            annotation_text="piQTL",
            annotation_position="right",
            row=1,
            col=1,
        )

    # Add reference markers for eQTL subplot
    if f"eQTL_global" in markers:
        fig.add_hline(
            y=markers["eQTL_global"][metric_col],
            line_dash="dash",
            line_color="red",
            line_width=2,
            annotation_text="Total eQTL",
            annotation_position="right",
            row=1,
            col=2,
        )
    if f"eQTL_piQTL_target_subset" in markers:
        fig.add_hline(
            y=markers["eQTL_piQTL_target_subset"][metric_col],
            line_dash="dot",
            line_color="green",
            line_width=2,
            annotation_text="44 piQTL-restricted eQTL",
            annotation_position="right",
            row=1,
            col=2,
        )
    if f"piQTL_global" in markers:
        fig.add_hline(
            y=markers["piQTL_global"][metric_col],
            line_dash="solid",
            line_color="purple",
            line_width=2,
            annotation_text="piQTL",
            annotation_position="right",
            row=1,
            col=2,
        )

    fig.update_yaxes(title_text=y_title, row=1, col=1)
    fig.update_yaxes(title_text=y_title, row=1, col=2)

    fig.update_layout(
        title_text=panel_title,
        title_font_size=18,
        height=600,
        width=1200,
        hovermode="closest",
    )

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Generate PNG and SVG only
    fig.write_image(str(output_path.with_suffix(".png")))
    fig.write_image(str(output_path.with_suffix(".svg")))

    print(f"Saved: {output_path.with_suffix('.png')}")
    print(f"Saved: {output_path.with_suffix('.svg')}")


def main():
    args = parse_args()

    print("Loading simulation data...")
    pqtl_sim, eqtl_sim = load_simulation_data(args.pqtl_sim, args.eqtl_sim)

    print("Loading reference markers...")
    markers = load_reference_markers(args.summary_table)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 72)
    print("Creating Figure Panels")
    print("=" * 72)

    # Panel 1: Mean difference
    print("\n[1/4] Mean difference (cis - trans)...")
    create_panel_figure(
        pqtl_sim=pqtl_sim,
        eqtl_sim=eqtl_sim,
        markers=markers,
        metric_col="mean_diff",
        y_title="Mean Difference (cis - trans)",
        panel_title="Panel 1: Mean Effect-Size Difference (cis - trans)",
        output_path=output_dir / "panel1_mean_diff",
    )

    # Panel 2: Median difference
    print("\n[2/4] Median difference (cis - trans)...")
    create_panel_figure(
        pqtl_sim=pqtl_sim,
        eqtl_sim=eqtl_sim,
        markers=markers,
        metric_col="median_diff",
        y_title="Median Difference (cis - trans)",
        panel_title="Panel 2: Median Effect-Size Difference (cis - trans)",
        output_path=output_dir / "panel2_median_diff",
    )

    # Panel 3: cis/trans ratio
    print("\n[3/4] cis/trans ratio...")
    create_panel_figure(
        pqtl_sim=pqtl_sim,
        eqtl_sim=eqtl_sim,
        markers=markers,
        metric_col="cis_trans_ratio",
        y_title="cis/trans Ratio",
        panel_title="Panel 3: cis/trans Ratio Distribution",
        output_path=output_dir / "panel3_cis_trans_ratio",
    )

    # Panel 4: trans percentage
    print("\n[4/4] trans-QTL percentage...")
    create_panel_figure(
        pqtl_sim=pqtl_sim,
        eqtl_sim=eqtl_sim,
        markers=markers,
        metric_col="trans_pct",
        y_title="trans-QTL Percentage (%)",
        panel_title="Panel 4: trans-QTL Percentage Distribution",
        output_path=output_dir / "panel4_trans_pct",
    )

    print("\n" + "=" * 72)
    print("All panels created successfully!")
    print("=" * 72)
    print(f"\nOutput directory: {output_dir}")
    print("\nGenerated files:")
    print("  - panel1_mean_diff.{png,svg}")
    print("  - panel2_median_diff.{png,svg}")
    print("  - panel3_cis_trans_ratio.{png,svg}")
    print("  - panel4_trans_pct.{png,svg}")


if __name__ == "__main__":
    main()
