"""Plot boxplot for comparing cis and trans effect"""

import argparse
import os

import pandas as pd
import plotly.express as px


# Load the command line arguments
def load_args():
    parser = argparse.ArgumentParser(description="Plot boxplot for comparing cis")
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to the input file",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to the output file",
    )
    parser.add_argument(
        "--title",
        type=str,
        required=True,
        help="Title of the plot",
    )
    parser.add_argument(
        "--y_label",
        type=str,
        required=True,
        help="Y-axis label",
    )
    parser.add_argument(
        "--group_by",
        type=str,
        required=True,
        help="Column to group by",
    )
    parser.add_argument(
        "--y_value",
        type=str,
        required=True,
        help="Column to use for y-axis values",
    )
    return parser.parse_args()


# Load the input data
def load_data(input_file):
    data = pd.read_csv(input_file, sep="\t", header=0)
    return data


# Plot boxplot
def plot_boxplot(data, title, y_label, group_by, y_value):
    # Specify the x-axis order: cis then trans
    x_order = ["cis", "trans"]
    x_label = "cis_trans"

    # Calculate sample counts for each group
    sample_counts = data[group_by].value_counts()
    sample_counts_text = "<br>".join(
        [f"{group}: {count}" for group, count in sample_counts.items()]
    )

    # Create the boxplot
    fig = px.box(
        data,
        x=group_by,
        y=y_value,
        title=title,
        labels={y_value: y_label, group_by: x_label},
        category_orders={group_by: x_order},
    )

    # Specify the plot size
    fig.update_layout(
        width=400,
        height=600,
    )
    # Add sample counts as an annotation
    fig.add_annotation(
        text=f"Sample counts:<br>{sample_counts_text}",
        xref="paper",
        yref="paper",
        showarrow=False,
        align="left",
        font=dict(size=12),
    )
    return fig


# main function
def main():
    # Load the command line arguments
    args = load_args()

    # Load the input data
    data = load_data(args.input)

    # Plot boxplot
    fig = plot_boxplot(
        data,
        args.title,
        args.y_label,
        args.group_by,
        args.y_value,
    )

    # Save the figure
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    fig.write_image(args.output)
    ## Save the figure as svg
    fig.write_image(args.output.replace(".png", ".svg"))


if __name__ == "__main__":
    main()
