#!/usr/bin/env python3
"""
Random sampling simulation for analyzing cis and trans pQTL effect distributions.

This script runs multiple simulations where random proteins are sampled from the
full protein list, and cis vs trans pQTL effect sizes are compared.

For each simulation:
- Randomly sample N proteins from the 1,226-protein universe
- Extract pQTLs affecting those proteins
- Calculate effect size metrics (mean, median) and counts for cis vs trans pQTLs
- Compute differences and ratios

Output:
- Summary table with metrics for each simulation
- List of sampled proteins for each simulation
- 5-panel boxplot visualization
"""

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

warnings.filterwarnings("ignore")


def load_data(protein_file, pqtl_file):
    """
    Load protein list and pQTL data.

    Parameters
    ----------
    protein_file : str
        Path to protein CSV file (should have 'ORF' or 'protein' column)
    pqtl_file : str
        Path to pQTL TSV file (should have 'protein', 'abs_beta', 'cis_trans' columns)

    Returns
    -------
    proteins : list
        List of protein ORF identifiers
    pqtl_df : pd.DataFrame
        pQTL data with columns: protein, abs_beta, cis_trans
    """
    # Load protein list
    protein_df = pd.read_csv(protein_file)

    # Get the protein column (try 'protein', 'ORF', or 'protein' variant)
    if "protein" in protein_df.columns:
        protein_col = "protein"
    elif "ORF" in protein_df.columns:
        protein_col = "ORF"
    else:
        # Use the first column
        protein_col = protein_df.columns[0]

    proteins = protein_df[protein_col].unique().tolist()

    # Load pQTL data
    pqtl_df = pd.read_csv(pqtl_file, sep="\t")

    # Ensure required columns exist
    required_cols = ["protein", "abs_beta", "cis_trans"]
    for col in required_cols:
        if col not in pqtl_df.columns:
            raise ValueError(f"Missing required column '{col}' in pQTL file")

    return proteins, pqtl_df


def run_single_simulation(proteins, pqtl_df, n_proteins, random_state=None):
    """
    Run a single simulation: sample proteins and calculate metrics.

    Parameters
    ----------
    proteins : list
        List of all available proteins
    pqtl_df : pd.DataFrame
        pQTL data
    n_proteins : int
        Number of proteins to sample
    random_state : int, optional
        Random seed for reproducibility

    Returns
    -------
    dict
        Dictionary containing:
        - sampled_proteins : list of protein ORFs
        - mean_beta_cis : mean effect size for cis-pQTLs
        - median_beta_cis : median effect size for cis-pQTLs
        - mean_beta_trans : mean effect size for trans-pQTLs
        - median_beta_trans : median effect size for trans-pQTLs
        - count_cis : number of cis-pQTLs
        - count_trans : number of trans-pQTLs
        - mean_diff : difference in means (cis - trans)
        - median_diff : difference in medians (cis - trans)
        - cis_trans_ratio : ratio of cis to trans pQTL counts
        - cis_effect_sizes : raw cis abs_beta values from this simulation
        - trans_effect_sizes : raw trans abs_beta values from this simulation
    """
    rng = np.random.RandomState(random_state)

    # Sample proteins
    sampled_proteins = rng.choice(proteins, size=n_proteins, replace=False).tolist()

    # Filter pQTLs for sampled proteins
    filtered_pqtls = pqtl_df[pqtl_df["protein"].isin(sampled_proteins)].copy()

    # Separate cis and trans pQTLs
    cis_pqtls = filtered_pqtls[filtered_pqtls["cis_trans"] == "cis"]["abs_beta"]
    trans_pqtls = filtered_pqtls[filtered_pqtls["cis_trans"] == "trans"]["abs_beta"]

    # Calculate metrics
    metrics = {
        "sampled_proteins": sampled_proteins,
        "cis_effect_sizes": cis_pqtls.tolist(),
        "trans_effect_sizes": trans_pqtls.tolist(),
        "mean_beta_cis": cis_pqtls.mean() if len(cis_pqtls) > 0 else np.nan,
        "median_beta_cis": cis_pqtls.median() if len(cis_pqtls) > 0 else np.nan,
        "mean_beta_trans": trans_pqtls.mean() if len(trans_pqtls) > 0 else np.nan,
        "median_beta_trans": trans_pqtls.median() if len(trans_pqtls) > 0 else np.nan,
        "count_cis": len(cis_pqtls),
        "count_trans": len(trans_pqtls),
        "mean_diff": (
            (cis_pqtls.mean() - trans_pqtls.mean())
            if len(cis_pqtls) > 0 and len(trans_pqtls) > 0
            else np.nan
        ),
        "median_diff": (
            (cis_pqtls.median() - trans_pqtls.median())
            if len(cis_pqtls) > 0 and len(trans_pqtls) > 0
            else np.nan
        ),
    }

    # Calculate ratio (handle division by zero)
    if metrics["count_trans"] > 0:
        metrics["cis_trans_ratio"] = metrics["count_cis"] / metrics["count_trans"]
    else:
        metrics["cis_trans_ratio"] = np.nan

    return metrics


def create_checkpoint_effect_size_boxplot(
    cis_effect_sizes,
    trans_effect_sizes,
    output_dir,
    checkpoint_num,
    n_proteins,
):
    """
    Create a single-panel grouped boxplot for one checkpoint simulation.

    Parameters
    ----------
    cis_effect_sizes : list[float]
        cis-pQTL abs_beta values from one simulation
    trans_effect_sizes : list[float]
        trans-pQTL abs_beta values from one simulation
    output_dir : str or Path
        Base output directory
    checkpoint_num : int
        Simulation index checkpoint (e.g., 100, 200, ...)
    n_proteins : int
        Number of proteins sampled in the simulation
    """
    checkpoints_dir = Path(output_dir) / "checkpoints"
    checkpoints_dir.mkdir(parents=True, exist_ok=True)

    fig = go.Figure()

    fig.add_trace(
        go.Box(
            y=cis_effect_sizes,
            name="cis",
            marker_color="#1f77b4",
            boxmean=False,
        )
    )

    fig.add_trace(
        go.Box(
            y=trans_effect_sizes,
            name="trans",
            marker_color="#ff7f0e",
            boxmean=False,
        )
    )

    fig.update_layout(
        title_text=(
            "Checkpoint Effect Size Distribution "
            f"(simulation {checkpoint_num}, {n_proteins} random proteins)"
        ),
        title_font_size=16,
        yaxis_title="abs_beta",
        xaxis_title="pQTL Type",
        boxmode="group",
        showlegend=False,
        height=500,
        width=400,
        hovermode="closest",
    )

    png_file = checkpoints_dir / f"boxplot_figure_cp{checkpoint_num}.png"
    fig.write_image(str(png_file))

    svg_file = checkpoints_dir / f"boxplot_figure_cp{checkpoint_num}.svg"
    fig.write_image(str(svg_file))

    print(f"  Saved checkpoint plot: {png_file}")
    print(f"  Saved checkpoint plot: {svg_file}")


def run_simulation(
    proteins,
    pqtl_df,
    n_proteins,
    n_simulations,
    random_seed=None,
    output_dir=None,
    checkpoint_interval=100,
):
    """
    Run multiple simulations.

    Parameters
    ----------
    proteins : list
        List of all available proteins
    pqtl_df : pd.DataFrame
        pQTL data
    n_proteins : int
        Number of proteins to sample per simulation
    n_simulations : int
        Number of simulations to run
    random_seed : int, optional
        Random seed for reproducibility
    output_dir : str or Path, optional
        Output directory where checkpoint figures will be saved
    checkpoint_interval : int, optional
        Save checkpoint figures every N simulations (default: 100)

    Returns
    -------
    summary_df : pd.DataFrame
        Summary metrics for each simulation
    protein_samples_df : pd.DataFrame
        Sampled protein names for each simulation
    """
    # Run simulations
    all_results = []
    protein_samples_list = []

    print(
        f"Running {n_simulations} simulations, sampling {n_proteins} proteins each..."
    )

    for i in range(n_simulations):
        if (i + 1) % 100 == 0:
            print(f"  Completed {i + 1}/{n_simulations} simulations")

        # Use different random state for each simulation
        sim_seed = None if random_seed is None else random_seed + i
        result = run_single_simulation(proteins, pqtl_df, n_proteins, sim_seed)

        sampled_proteins = result.pop("sampled_proteins")
        cis_effect_sizes = result.pop("cis_effect_sizes")
        trans_effect_sizes = result.pop("trans_effect_sizes")

        all_results.append(result)
        protein_samples_list.append(sampled_proteins)

        if (
            output_dir is not None
            and checkpoint_interval > 0
            and (i + 1) % checkpoint_interval == 0
        ):
            create_checkpoint_effect_size_boxplot(
                cis_effect_sizes=cis_effect_sizes,
                trans_effect_sizes=trans_effect_sizes,
                output_dir=output_dir,
                checkpoint_num=i + 1,
                n_proteins=n_proteins,
            )

    print(f"  Completed {n_simulations}/{n_simulations} simulations")

    # Create summary dataframe
    summary_df = pd.DataFrame(all_results)

    # Create protein samples dataframe (padded to match simulation count)
    # Each row is a simulation, each column is a protein
    max_proteins = max(len(p) for p in protein_samples_list)
    protein_samples_data = []
    for proteins_in_sim in protein_samples_list:
        # Pad with empty strings if needed
        padded = proteins_in_sim + [""] * (max_proteins - len(proteins_in_sim))
        protein_samples_data.append(padded)

    protein_samples_df = pd.DataFrame(
        protein_samples_data, columns=[f"protein_{i+1}" for i in range(max_proteins)]
    )

    return summary_df, protein_samples_df


def create_boxplot_figure(summary_df, output_dir, n_proteins=44, n_simulations=1000):
    """
    Create a 5-panel boxplot figure.

    Parameters
    ----------
    summary_df : pd.DataFrame
        Summary metrics from simulations
    output_dir : str
        Output directory for figure files
    n_proteins : int, optional
        Number of proteins sampled per simulation (default: 44)
    n_simulations : int, optional
        Number of simulations run (default: 1000)
    """
    fig = make_subplots(
        rows=1,
        cols=5,
        subplot_titles=(
            "Mean Effect Size Diff<br>(cis - trans)",
            "Median Effect Size Diff<br>(cis - trans)",
            "cis-pQTL Count",
            "trans-pQTL Count",
            "cis/trans Ratio",
        ),
        horizontal_spacing=0.08,
    )

    # Data for each subplot
    data_dict = {
        "mean_diff": summary_df["mean_diff"].dropna(),
        "median_diff": summary_df["median_diff"].dropna(),
        "count_cis": summary_df["count_cis"],
        "count_trans": summary_df["count_trans"],
        "cis_trans_ratio": summary_df["cis_trans_ratio"].dropna(),
    }

    y_axes = ["y", "y2", "y3", "y4", "y5"]
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]

    # Add boxplots
    for col_idx, (key, data) in enumerate(data_dict.items(), 1):
        fig.add_trace(
            go.Box(
                y=data,
                name=key,
                marker_color=colors[col_idx - 1],
                # boxmean="sd",
                yaxis=y_axes[col_idx - 1],
            ),
            row=1,
            col=col_idx,
        )

    # Update layout
    title_text = f"pQTL Simulation Results: Cis vs Trans Effect Size Metrics ({n_proteins} random proteins, {n_simulations} iterations)"
    fig.update_layout(
        height=500,
        width=1400,
        showlegend=False,
        title_text=title_text,
        title_font_size=16,
        hovermode="closest",
    )

    # Update y-axes
    y_axis_titles = ["Mean Diff", "Median Diff", "Count", "Count", "Ratio"]
    for col_idx, title in enumerate(y_axis_titles, 1):
        fig.update_yaxes(title_text=title, row=1, col=col_idx)

    # Save figures
    output_path = Path(output_dir)

    # HTML (interactive)
    html_file = output_path / "boxplot_figure.html"
    fig.write_html(str(html_file))
    print(f"  Saved interactive plot: {html_file}")

    # PNG (static)
    png_file = output_path / "boxplot_figure.png"
    fig.write_image(str(png_file), width=1400, height=500)
    print(f"  Saved PNG plot: {png_file}")

    # SVG (static)
    svg_file = output_path / "boxplot_figure.svg"
    fig.write_image(str(svg_file), width=1400, height=500)
    print(f"  Saved SVG plot: {svg_file}")


def save_summary_statistics(
    summary_df, output_dir, n_proteins, n_simulations, random_seed
):
    """
    Save summary statistics to a text file.

    Parameters
    ----------
    summary_df : pd.DataFrame
        Summary metrics from simulations
    output_dir : str
        Output directory
    n_proteins : int
        Number of proteins sampled per simulation
    n_simulations : int
        Number of simulations run
    random_seed : int or None
        Random seed used
    """
    output_path = Path(output_dir)
    stats_file = output_path / "summary_statistics.txt"

    with open(stats_file, "w") as f:
        f.write("=" * 70 + "\n")
        f.write("RANDOM SAMPLING pQTL SIMULATION - SUMMARY STATISTICS\n")
        f.write("=" * 70 + "\n\n")

        f.write("SIMULATION PARAMETERS:\n")
        f.write(f"  Number of simulations: {n_simulations}\n")
        f.write(f"  Proteins sampled per simulation: {n_proteins}\n")
        f.write(f"  Total proteins in universe: {1225}\n")  # Excluding header
        f.write(f"  Total pQTLs available: {len(summary_df)}\n")  # This will be updated
        f.write(
            f"  Random seed: {random_seed if random_seed is not None else 'Not set (non-reproducible)'}\n"
        )
        f.write("\n")

        metrics = [
            ("mean_diff", "Mean Effect Size Difference (cis - trans)"),
            ("median_diff", "Median Effect Size Difference (cis - trans)"),
            ("count_cis", "cis-pQTL Count"),
            ("count_trans", "trans-pQTL Count"),
            ("cis_trans_ratio", "cis/trans pQTL Ratio"),
        ]

        f.write("SUMMARY STATISTICS FOR KEY METRICS:\n")
        f.write("-" * 70 + "\n")

        for col, label in metrics:
            data = summary_df[col].dropna()
            f.write(f"\n{label}:\n")
            f.write(f"  N (simulations with data):     {len(data)}\n")
            f.write(f"  Mean:                          {data.mean():.6f}\n")
            f.write(f"  Median:                        {data.median():.6f}\n")
            f.write(f"  Std Dev:                       {data.std():.6f}\n")
            f.write(f"  Min:                           {data.min():.6f}\n")
            f.write(f"  Max:                           {data.max():.6f}\n")
            f.write(f"  Q1 (25th percentile):          {data.quantile(0.25):.6f}\n")
            f.write(f"  Q3 (75th percentile):          {data.quantile(0.75):.6f}\n")

        f.write("\n" + "=" * 70 + "\n")

    print(f"  Saved summary statistics: {stats_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Random sampling simulation for cis/trans pQTL effect analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Default run: 1,000 simulations sampling 44 proteins each
  python random_sampling_pQTL_simulation.py

  # Custom parameters with reproducible seed
  python random_sampling_pQTL_simulation.py \\
    --n-proteins 50 \\
    --n-simulations 500 \\
    --random-seed 42

  # Specify data and output directories
  python random_sampling_pQTL_simulation.py \\
    --protein-file /path/to/proteins.csv \\
    --pqtl-file /path/to/pqtls.tsv \\
    --output-dir /path/to/output
        """,
    )

    parser.add_argument(
        "--protein-file",
        default="data/Jakobson_tested_genes.csv",
        help="Path to protein list CSV file (default: data/Jakobson_tested_genes.csv)",
    )

    parser.add_argument(
        "--pqtl-file",
        default="data/Jakobson_pQTLs.tsv",
        help="Path to pQTL TSV file (default: data/Jakobson_pQTLs.tsv)",
    )

    parser.add_argument(
        "--n-proteins",
        type=int,
        default=44,
        help="Number of proteins to sample per simulation (default: 44)",
    )

    parser.add_argument(
        "--n-simulations",
        type=int,
        default=1000,
        help="Number of simulations to run (default: 1000)",
    )

    parser.add_argument(
        "--random-seed",
        type=int,
        default=None,
        help="Random seed for reproducibility (default: None, non-reproducible)",
    )

    parser.add_argument(
        "--output-dir",
        default="simulation_results",
        help="Output directory for results (default: simulation_results)",
    )

    args = parser.parse_args()

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {output_dir}")

    # Load data
    print(f"\nLoading data...")
    print(f"  Protein file: {args.protein_file}")
    print(f"  pQTL file: {args.pqtl_file}")
    proteins, pqtl_df = load_data(args.protein_file, args.pqtl_file)
    print(f"  Loaded {len(proteins)} proteins")
    print(f"  Loaded {len(pqtl_df)} pQTLs")

    # Run simulations
    print(f"\nRunning simulations...")
    summary_df, protein_samples_df = run_simulation(
        proteins,
        pqtl_df,
        args.n_proteins,
        args.n_simulations,
        args.random_seed,
        output_dir=output_dir,
        checkpoint_interval=100,
    )

    # Save results
    print(f"\nSaving results...")

    # Save summary table
    summary_file = output_dir / "simulation_summary.tsv"
    summary_df.to_csv(summary_file, sep="\t", index=False)
    print(f"  Saved summary table: {summary_file}")
    print(f"    Shape: {summary_df.shape[0]} rows × {summary_df.shape[1]} columns")

    # Save protein samples table
    protein_samples_file = output_dir / "protein_samples.tsv"
    protein_samples_df.to_csv(protein_samples_file, sep="\t", index=False)
    print(f"  Saved protein samples table: {protein_samples_file}")
    print(
        f"    Shape: {protein_samples_df.shape[0]} rows × {protein_samples_df.shape[1]} columns"
    )

    # Create visualization
    print(f"\nCreating figures...")
    create_boxplot_figure(summary_df, output_dir, args.n_proteins, args.n_simulations)

    # Save summary statistics
    print(f"\nSaving summary statistics...")
    save_summary_statistics(
        summary_df, output_dir, args.n_proteins, args.n_simulations, args.random_seed
    )

    print(f"\n" + "=" * 70)
    print(f"Simulation completed successfully!")
    print(f"Results saved to: {output_dir}")
    print(f"=" * 70)


if __name__ == "__main__":
    main()
