"""t-test for cis and trans QTLs"""

import argparse
import os

import pandas as pd
from scipy.stats import ttest_ind


# Load the command line arguments
def load_args():
    parser = argparse.ArgumentParser(description="Perform t-test for QTLs")

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
    parser.add_argument(
        "--title",
        type=str,
        required=True,
        help="Title of the results",
    )
    return parser.parse_args()


# Load the input data
def load_data(input_file):
    data = pd.read_csv(input_file, sep="\t", header=0)
    return data


# Perform t-test
def perform_ttest(data, group_by, y_value):
    # Split the data into cis and trans groups
    cis_group = data[data[group_by] == "cis"][y_value]
    trans_group = data[data[group_by] == "trans"][y_value]
    # Perform t-test
    t_stat, p_value = ttest_ind(cis_group, trans_group, equal_var=False)
    return t_stat, p_value


# Save the results
def save_results(output_file, t_stat, p_value, title):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, "w") as f:
        f.write(f"{title}\n")
        f.write(f"t-statistic: {t_stat}\n")
        f.write(f"p-value: {p_value}\n")


# Main function
def main():
    # Load the command line arguments
    args = load_args()

    # Load the input data
    data = load_data(args.input)

    # Perform t-test
    t_stat, p_value = perform_ttest(data, args.group_by, args.y_value)

    # Save the results
    save_results(args.output, t_stat, p_value, args.title)


if __name__ == "__main__":
    main()
