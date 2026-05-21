#!/usr/bin/env python3
"""
Random sampling simulation for analyzing cis and trans eQTL effect distributions.

This script runs multiple simulations where random genes are sampled from the
full gene list, and cis vs trans eQTL effect sizes are compared.

For each simulation:
- Randomly sample N genes from the unique gene universe
- Extract eQTLs affecting those genes
- Calculate effect size metrics (mean, median) and counts for cis vs trans eQTLs
- Compute differences and ratios

Output:
- Summary table with metrics for each simulation
- List of sampled genes for each simulation
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


def load_data(eqtl_file):
    """
    Load eQTL data and extract unique gene list.

    Parameters
    ----------
    eqtl_file : str
        Path to eQTL TSV file (should have 'gene', 'beta', 'cis_trans' columns)

    Returns
    -------
    genes : list
        List of unique gene identifiers
    eqtl_df : pd.DataFrame
        eQTL data with columns: gene, abs_beta, cis_trans
    """
    # Load eQTL data
    eqtl_df = pd.read_csv(eqtl_file, sep="\t")

    # Ensure required columns exist
    required_cols = ["gene", "beta", "cis_trans"]
    for col in required_cols:
        if col not in eqtl_df.columns:
            raise ValueError(f"Missing required column '{col}' in eQTL file")

    # Compute abs_beta from beta
    eqtl_df["abs_beta"] = eqtl_df["beta"].abs()

    # Extract unique genes
    genes = eqtl_df["gene"].unique().tolist()

    return genes, eqtl_df


def run_single_simulation(genes, eqtl_df, n_genes, random_state=None):
    """
    Run a single simulation: sample genes and calculate metrics.

    Parameters
    ----------
    genes : list
        List of all available genes
    eqtl_df : pd.DataFrame
        eQTL data
    n_genes : int
        Number of genes to sample
    random_state : int, optional
        Random seed for reproducibility

    Returns
    -------
    dict
        Dictionary containing:
        - sampled_genes : list of gene identifiers
        - mean_beta_cis : mean effect size for cis-eQTLs
        - median_beta_cis : median effect size for cis-eQTLs
        - mean_beta_trans : mean effect size for trans-eQTLs
        - median_beta_trans : median effect size for trans-eQTLs
        - count_cis : number of cis-eQTLs
        - count_trans : number of trans-eQTLs
        - mean_diff : difference in means (cis - trans)
        - median_diff : difference in medians (cis - trans)
        - cis_trans_ratio : ratio of cis to trans eQTL counts
        - cis_effect_sizes : raw cis abs_beta values from this simulation
        - trans_effect_sizes : raw trans abs_beta values from this simulation
    """
    rng = np.random.RandomState(random_state)

    # Sample genes
    sampled_genes = rng.choice(genes, size=n_genes, replace=False).tolist()

    # Filter eQTLs for sampled genes
    filtered_eqtls = eqtl_df[eqtl_df["gene"].isin(sampled_genes)].copy()

    # Separate cis and trans eQTLs
    cis_eqtls = filtered_eqtls[filtered_eqtls["cis_trans"] == "cis"]["abs_beta"]
    trans_eqtls = filtered_eqtls[filtered_eqtls["cis_trans"] == "trans"]["abs_beta"]

    # Calculate metrics
    metrics = {
        "sampled_genes": sampled_genes,
        "cis_effect_sizes": cis_eqtls.tolist(),
        "trans_effect_sizes": trans_eqtls.tolist(),
        "mean_beta_cis": cis_eqtls.mean() if len(cis_eqtls) > 0 else np.nan,
        "median_beta_cis": cis_eqtls.median() if len(cis_eqtls) > 0 else np.nan,
        "mean_beta_trans": trans_eqtls.mean() if len(trans_eqtls) > 0 else np.nan,
        "median_beta_trans": trans_eqtls.median() if len(trans_eqtls) > 0 else np.nan,
        "count_cis": len(cis_eqtls),
        "count_trans": len(trans_eqtls),
        "mean_diff": (
            (cis_eqtls.mean() - trans_eqtls.mean())
            if len(cis_eqtls) > 0 and len(trans_eqtls) > 0
            else np.nan
        ),
        "median_diff": (
            (cis_eqtls.median() - trans_eqtls.median())
            if len(cis_eqtls) > 0 and len(trans_eqtls) > 0
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
    n_genes,
):
    """
    Create a single-panel grouped boxplot for one checkpoint simulation.

    Parameters
    ----------
    cis_effect_sizes : list[float]
        cis-eQTL abs_beta values from one simulation
    trans_effect_sizes : list[float]
        trans-eQTL abs_beta values from one simulation
    output_dir : str or Path
        Base output directory
    checkpoint_num : int
        Simulation index checkpoint (e.g., 100, 200, ...)
    n_genes : int
        Number of genes sampled in the simulation
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
            f"(simulation {checkpoint_num}, {n_genes} random genes)"
        ),
        title_font_size=16,
        yaxis_title="abs_beta",
        xaxis_title="eQTL Type",
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
    genes,
    eqtl_df,
    n_genes,
    n_simulations,
    random_seed=None,
    output_dir=None,
    checkpoint_interval=100,
):
    """
    Run multiple simulations.

    Parameters
    ----------
    genes : list
        List of all available genes
    eqtl_df : pd.DataFrame
        eQTL data
    n_genes : int
        Number of genes to sample per simulation
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
    gene_samples_df : pd.DataFrame
        Sampled gene names for each simulation
    """
    # Run simulations
    all_results = []
    gene_samples_list = []

    print(f"Running {n_simulations} simulations, sampling {n_genes} genes each...")

    for i in range(n_simulations):
        if (i + 1) % 100 == 0:
            print(f"  Completed {i + 1}/{n_simulations} simulations")

        # Use different random state for each simulation
        sim_seed = None if random_seed is None else random_seed + i
        result = run_single_simulation(genes, eqtl_df, n_genes, sim_seed)

        sampled_genes = result.pop("sampled_genes")
        cis_effect_sizes = result.pop("cis_effect_sizes")
        trans_effect_sizes = result.pop("trans_effect_sizes")

        all_results.append(result)
        gene_samples_list.append(sampled_genes)

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
                n_genes=n_genes,
            )

    print(f"  Completed {n_simulations}/{n_simulations} simulations")

    # Create summary dataframe
    summary_df = pd.DataFrame(all_results)

    # Create gene samples dataframe (padded to match simulation count)
    # Each row is a simulation, each column is a gene
    max_genes = max(len(g) for g in gene_samples_list)
    gene_samples_data = []
    for genes_in_sim in gene_samples_list:
        # Pad with empty strings if needed
        padded = genes_in_sim + [""] * (max_genes - len(genes_in_sim))
        gene_samples_data.append(padded)

    gene_samples_df = pd.DataFrame(
        gene_samples_data, columns=[f"gene_{i+1}" for i in range(max_genes)]
    )

    return summary_df, gene_samples_df


def create_boxplot_figure(summary_df, output_dir, n_genes=44, n_simulations=1000):
    """
    Create a 5-panel boxplot figure.

    Parameters
    ----------
    summary_df : pd.DataFrame
        Summary metrics from simulations
    output_dir : str
        Output directory for figure files
    n_genes : int, optional
        Number of genes sampled per simulation (default: 44)
    n_simulations : int, optional
        Number of simulations run (default: 1000)
    """
    fig = make_subplots(
        rows=1,
        cols=5,
        subplot_titles=(
            "Mean Effect Size Diff<br>(cis - trans)",
            "Median Effect Size Diff<br>(cis - trans)",
            "cis-eQTL Count",
            "trans-eQTL Count",
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
    title_text = f"eQTL Simulation Results: Cis vs Trans Effect Size Metrics ({n_genes} random genes, {n_simulations} iterations)"
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
    summary_df,
    output_dir,
    n_genes,
    n_simulations,
    random_seed,
    total_genes,
    total_eqtls,
):
    """
    Save summary statistics to a text file.

    Parameters
    ----------
    summary_df : pd.DataFrame
        Summary metrics from simulations
    output_dir : str
        Output directory
    n_genes : int
        Number of genes sampled per simulation
    n_simulations : int
        Number of simulations run
    random_seed : int or None
        Random seed used
    total_genes : int
        Total unique genes in the universe
    total_eqtls : int
        Total eQTLs available
    """
    output_path = Path(output_dir)
    stats_file = output_path / "summary_statistics.txt"

    with open(stats_file, "w") as f:
        f.write("=" * 70 + "\n")
        f.write("RANDOM SAMPLING eQTL SIMULATION - SUMMARY STATISTICS\n")
        f.write("=" * 70 + "\n\n")

        f.write("SIMULATION PARAMETERS:\n")
        f.write(f"  Number of simulations: {n_simulations}\n")
        f.write(f"  Genes sampled per simulation: {n_genes}\n")
        f.write(f"  Total genes in universe: {total_genes}\n")
        f.write(f"  Total eQTLs available: {total_eqtls}\n")
        f.write(
            f"  Random seed: {random_seed if random_seed is not None else 'Not set (non-reproducible)'}\n"
        )
        f.write("\n")

        metrics = [
            ("mean_diff", "Mean Effect Size Difference (cis - trans)"),
            ("median_diff", "Median Effect Size Difference (cis - trans)"),
            ("count_cis", "cis-eQTL Count"),
            ("count_trans", "trans-eQTL Count"),
            ("cis_trans_ratio", "cis/trans eQTL Ratio"),
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
        description="Random sampling simulation for cis/trans eQTL effect analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Default run: 1,000 simulations sampling 44 genes each
  python random_sampling_eQTL_simulation.py

  # Custom parameters with reproducible seed
  python random_sampling_eQTL_simulation.py \\
    --n-genes 50 \\
    --n-simulations 500 \\
    --random-seed 42

  # Specify data and output directories
  python random_sampling_eQTL_simulation.py \\
    --eqtl-file /path/to/eqtls.tsv \\
    --output-dir /path/to/output
        """,
    )

    parser.add_argument(
        "--eqtl-file",
        default="data/Albert_eQTLs.tsv",
        help="Path to eQTL TSV file (default: data/Albert_eQTLs.tsv)",
    )

    parser.add_argument(
        "--n-genes",
        type=int,
        default=44,
        help="Number of genes to sample per simulation (default: 44)",
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
        default="out/simulation_results_eQTL",
        help="Output directory for results (default: out/simulation_results_eQTL)",
    )

    args = parser.parse_args()

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {output_dir}")

    # Load data
    print(f"\nLoading data...")
    print(f"  eQTL file: {args.eqtl_file}")
    genes, eqtl_df = load_data(args.eqtl_file)
    print(f"  Loaded {len(genes)} unique genes")
    print(f"  Loaded {len(eqtl_df)} eQTLs")

    # Run simulations
    print(f"\nRunning simulations...")
    summary_df, gene_samples_df = run_simulation(
        genes,
        eqtl_df,
        args.n_genes,
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

    # Save gene samples table
    gene_samples_file = output_dir / "gene_samples.tsv"
    gene_samples_df.to_csv(gene_samples_file, sep="\t", index=False)
    print(f"  Saved gene samples table: {gene_samples_file}")
    print(
        f"    Shape: {gene_samples_df.shape[0]} rows × {gene_samples_df.shape[1]} columns"
    )

    # Create visualization
    print(f"\nCreating figures...")
    create_boxplot_figure(summary_df, output_dir, args.n_genes, args.n_simulations)

    # Save summary statistics
    print(f"\nSaving summary statistics...")
    save_summary_statistics(
        summary_df,
        output_dir,
        args.n_genes,
        args.n_simulations,
        args.random_seed,
        len(genes),
        len(eqtl_df),
    )

    print(f"\n" + "=" * 70)
    print(f"Simulation completed successfully!")
    print(f"Results saved to: {output_dir}")
    print(f"=" * 70)


if __name__ == "__main__":
    main()
