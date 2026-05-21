"""Plot boxplot for Jakobson pQTLs filtered to piQTL genes and perform t-test"""

import argparse
import os

import pandas as pd
import plotly.express as px
from scipy.stats import ttest_ind


def load_args():
    """Load command line arguments"""
    parser = argparse.ArgumentParser(
        description="Plot boxplot for Jakobson pQTLs filtered to piQTL genes"
    )
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to the filtered pQTL input file (TSV format)",
    )
    parser.add_argument(
        "--output_plot",
        type=str,
        required=True,
        help="Path to the output plot file (PNG)",
    )
    parser.add_argument(
        "--output_stats",
        type=str,
        required=True,
        help="Path to the output statistics file (TXT)",
    )
    return parser.parse_args()


def load_data(input_file):
    """Load the filtered pQTL data"""
    data = pd.read_csv(input_file, sep="\t", header=0)
    print(f"Loaded {len(data)} records from {input_file}")
    return data


def perform_ttest(data):
    """Perform t-test comparing abs_beta between cis and trans"""
    cis_data = data[data["cis_trans"] == "cis"]["abs_beta"].values
    trans_data = data[data["cis_trans"] == "trans"]["abs_beta"].values

    # Perform independent samples t-test
    t_stat, p_value = ttest_ind(cis_data, trans_data)

    # Calculate summary statistics
    cis_mean = cis_data.mean()
    cis_std = cis_data.std()
    cis_count = len(cis_data)

    trans_mean = trans_data.mean()
    trans_std = trans_data.std()
    trans_count = len(trans_data)

    return {
        "cis": {
            "n": cis_count,
            "mean": cis_mean,
            "std": cis_std,
        },
        "trans": {
            "n": trans_count,
            "mean": trans_mean,
            "std": trans_std,
        },
        "t_stat": t_stat,
        "p_value": p_value,
    }


def save_statistics(stats, output_file):
    """Save t-test statistics to file"""
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with open(output_file, "w") as f:
        f.write("T-Test Results: Comparing cis vs trans pQTL effect sizes (abs_beta)\n")
        f.write("=" * 70 + "\n\n")

        f.write("CIS pQTLs:\n")
        f.write(f"  Sample size: {stats['cis']['n']}\n")
        f.write(f"  Mean effect size: {stats['cis']['mean']:.6f}\n")
        f.write(f"  Std deviation: {stats['cis']['std']:.6f}\n\n")

        f.write("TRANS pQTLs:\n")
        f.write(f"  Sample size: {stats['trans']['n']}\n")
        f.write(f"  Mean effect size: {stats['trans']['mean']:.6f}\n")
        f.write(f"  Std deviation: {stats['trans']['std']:.6f}\n\n")

        f.write("Independent samples t-test:\n")
        f.write(f"  t-statistic: {stats['t_stat']:.6f}\n")
        f.write(f"  p-value: {stats['p_value']:.6e}\n")

    print(f"Saved statistics to {output_file}")


def plot_boxplot(data, output_plot):
    """Create and save boxplot"""
    # Specify the x-axis order: cis then trans
    x_order = ["cis", "trans"]

    # Calculate sample counts for each group
    sample_counts = data["cis_trans"].value_counts().to_dict()
    sample_counts_text = "<br>".join(
        [f"{group}: {count}" for group, count in sorted(sample_counts.items())]
    )

    # Create the boxplot
    fig = px.box(
        data,
        x="cis_trans",
        y="abs_beta",
        title="Jakobson pQTL Effect Sizes (Filtered to piQTL Target Genes)",
        labels={"abs_beta": "Absolute Effect Size (|β|)", "cis_trans": "QTL Type"},
        category_orders={"cis_trans": x_order},
        points="all",  # Show all data points
        color="cis_trans",
    )

    # Specify the plot size
    fig.update_layout(
        width=600,
        height=700,
        showlegend=False,
    )

    # Add sample counts as an annotation
    fig.add_annotation(
        text=f"Sample counts:<br>{sample_counts_text}",
        xref="paper",
        yref="paper",
        x=0.80,
        y=0.98,
        showarrow=False,
        align="left",
        font=dict(size=11),
        bgcolor="rgba(255, 255, 255, 0.8)",
        bordercolor="black",
        borderwidth=1,
    )

    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_plot), exist_ok=True)

    # Save the figure as PNG
    fig.write_image(output_plot)
    print(f"Saved plot to {output_plot}")

    # Save the figure as SVG
    svg_file = output_plot.replace(".png", ".svg")
    fig.write_image(svg_file)
    print(f"Saved plot to {svg_file}")


def main():
    # Load arguments
    args = load_args()

    # Load data
    data = load_data(args.input)

    # Perform t-test
    stats = perform_ttest(data)

    # Save statistics
    save_statistics(stats, args.output_stats)

    # Create and save boxplot
    plot_boxplot(data, args.output_plot)

    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()
