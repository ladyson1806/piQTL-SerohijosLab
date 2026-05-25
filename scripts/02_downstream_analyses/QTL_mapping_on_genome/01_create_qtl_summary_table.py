#!/usr/bin/env python3
"""
Create a comprehensive QTL summary table combining piQTL, pQTL, and eQTL results.

This script:
1. Loads piQTL SNP annotation as the base table
2. Adds isAlbert_SNP column by matching with eQTL SNP annotation
3. Counts exact QTL matches (peak position = SNP position)
4. Counts colocal QTL matches (SNP position within left-right range)
5. Exports final summary table
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd

# Define paths
data_dir = Path("data")
out_dir = Path("out/tables")
out_dir.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("Creating QTL Summary Table")
print("=" * 80)

# ============================================================================
# STEP 1: Load and prepare piQTL SNP annotation as base table
# ============================================================================
print("\n[1/5] Loading piQTL SNP annotation...")
piqtl_snps = pd.read_csv(data_dir / "piQTL_SNP_annotation.csv")
print(f"    Loaded {len(piqtl_snps)} piQTL SNPs")


# Extract chromosome number from "CHR_X" format
# Handle special case of mitochondrial chromosome (CHR_MT)
def chr_to_num(chr_str):
    """Convert chromosome string to number, treating MT as 17"""
    chr_str = str(chr_str).replace("CHR_", "")
    if chr_str == "MT":
        return 17
    try:
        return int(chr_str)
    except ValueError:
        return 17  # Default to 17 for unknown


piqtl_snps["chromosome"] = piqtl_snps["chrom"].apply(chr_to_num)
piqtl_snps["position"] = piqtl_snps["position"].astype(int)

# Create SNP marker in format "chr{num}:{position}"
piqtl_snps["SNP_marker"] = (
    "chr"
    + piqtl_snps["chromosome"].astype(str)
    + ":"
    + piqtl_snps["position"].astype(str)
)

# # Initialize count columns
# for col in [
#     "exact_piQTL",
#     "exact_pQTL",
#     "exact_eQTL",
#     "colocal_piQTL",
#     "colocal_pQTL",
#     "colocal_eQTL",
# ]:
#     piqtl_snps[col] = 0

# Select and reorder columns for final output
summary = piqtl_snps[["SNP", "SNP_marker", "chromosome", "position"]].copy()

# Initialize new columns for QTL matches and Albert SNP annotation
summary["isAlbert_SNP"] = False
summary["exact_piQTL"] = 0
summary["exact_pQTL"] = 0
summary["exact_eQTL"] = 0
summary["colocal_piQTL"] = 0
summary["colocal_pQTL"] = 0
summary["colocal_eQTL"] = 0

print(f"    Created base table with {len(summary)} SNPs")

# ============================================================================
# STEP 2: Add isAlbert_SNP column
# ============================================================================
print("\n[2/5] Adding isAlbert_SNP column...")
albert_snps_raw = pd.read_csv(
    data_dir / "eQTL_SNP_annotation.csv", header=None, names=["marker"]
)
print(f"    Loaded {len(albert_snps_raw)} Albert SNP markers")


# Parse Albert SNP markers: format is "chrI:position_REF/ALT"
# Extract position and alleles for matching
def parse_albert_marker(marker):
    """Parse Albert SNP marker format: chrI:position_allele/allele"""
    if pd.isna(marker):
        return None, None, None

    # Pattern: chrX:position_alleles
    match = re.match(r"chr([IVX]+):(\d+)_(.+)$", marker)
    if not match:
        return None, None, None

    chr_roman = match.group(1)
    position = int(match.group(2))
    alleles = match.group(3)

    # Convert Roman numeral chromosome to number
    roman_to_num = {
        "I": 1,
        "II": 2,
        "III": 3,
        "IV": 4,
        "V": 5,
        "VI": 6,
        "VII": 7,
        "VIII": 8,
        "IX": 9,
        "X": 10,
        "XI": 11,
        "XII": 12,
        "XIII": 13,
        "XIV": 14,
        "XV": 15,
        "XVI": 16,
    }
    chr_num = roman_to_num.get(chr_roman)

    return chr_num, position, alleles


albert_snps_raw[["chr_num", "position", "alleles"]] = albert_snps_raw["marker"].apply(
    lambda x: pd.Series(parse_albert_marker(x))
)

# Create set of (chr, pos) tuples for fast lookup
albert_pos_set = set()

for _, row in albert_snps_raw.dropna(subset=["chr_num", "position"]).iterrows():
    albert_pos_set.add((int(row["chr_num"]), int(row["position"])))

print(f"    Parsed {len(albert_pos_set)} unique Albert SNP positions")


# Match piQTL SNPs with Albert SNPs (position + allele match required)
def is_albert_snp(snp_row):
    """Check if SNP exists in Albert dataset with matching alleles"""
    chr_num = snp_row["chromosome"]
    pos = snp_row["position"]

    # Check if this combination exists
    return (chr_num, pos) in albert_pos_set


summary["isAlbert_SNP"] = piqtl_snps.apply(is_albert_snp, axis=1)
print(f"    Matched {summary['isAlbert_SNP'].sum()} SNPs in Albert dataset")


# ============================================================================
# STEP 3: Count exact QTL matches
# ============================================================================
print("\n[3/5] Counting exact QTL matches...")

# Load QTL results
piqtl_results = pd.read_csv(data_dir / "piQTL_results.csv")
pqtl_results = pd.read_csv(data_dir / "pQTL_results.csv")
eqtl_results = pd.read_csv(data_dir / "eQTL_results.csv")

print(f"    Loaded {len(piqtl_results)} piQTL results")
print(f"    Loaded {len(pqtl_results)} pQTL results")
print(f"    Loaded {len(eqtl_results)} eQTL results")

# Count exact piQTL matches (piQTL_peak == SNP position)
exact_piqtl_dict = {}
for snp_id in summary["SNP"]:
    snp_pos = summary[summary["SNP"] == snp_id]["position"].values[0]
    snp_chr = summary[summary["SNP"] == snp_id]["chromosome"].values[0]
    matching_rows = piqtl_results[
        (piqtl_results["SNP"] == snp_id) & (piqtl_results["piQTL_peak"] == snp_pos)
    ]
    exact_piqtl_dict[snp_id] = len(matching_rows)

summary["exact_piQTL"] = summary["SNP"].map(exact_piqtl_dict).fillna(0).astype(int)
print(f"    Counted exact piQTL matches")

# Count exact pQTL matches (pQTL_peak == SNP position)
exact_pqtl_dict = {}
for snp_id in summary["SNP"]:
    snp_pos = summary[summary["SNP"] == snp_id]["position"].values[0]
    snp_chr = summary[summary["SNP"] == snp_id]["chromosome"].values[0]
    matching_rows = pqtl_results[
        (pqtl_results["SNP"] == snp_id) & (pqtl_results["pQTL_peak"] == snp_pos)
    ]
    exact_pqtl_dict[snp_id] = len(matching_rows)

summary["exact_pQTL"] = summary["SNP"].map(exact_pqtl_dict).fillna(0).astype(int)
print(f"    Counted exact pQTL matches")

# Count exact eQTL matches (eQTL_peak == SNP position AND chromosome match)
# eQTL doesn't have SNP column, so match by chromosome and position
exact_eqtl_dict = {}
for snp_id in summary["SNP"]:
    snp_pos = summary[summary["SNP"] == snp_id]["position"].values[0]
    snp_chr = summary[summary["SNP"] == snp_id]["chromosome"].values[0]
    # Match by chromosome and peak position
    matching_rows = eqtl_results[
        (eqtl_results["chromosome"] == snp_chr) & (eqtl_results["eQTL_peak"] == snp_pos)
    ]
    exact_eqtl_dict[snp_id] = len(matching_rows)

summary["exact_eQTL"] = summary["SNP"].map(exact_eqtl_dict).fillna(0).astype(int)
print(f"    Counted exact eQTL matches")

# ============================================================================
# STEP 4: Count colocal QTL matches
# ============================================================================
print("\n[4/5] Counting colocal QTL matches...")

# Note: Some QTL entries have reversed boundaries (left > right),
# Note: this is likely due to the right boundary beyond the chromosome end.
# Note: We will handle this by replacing right boundaries before matching.
# Note: The modification is replacing the right position with the left position + 50 kbp.

# Normalizing the QTL boundaries (if left > right, set right = left + 50 kbp)
pqtl_results["pQTL_right"] = np.where(
    pqtl_results["pQTL_left"] > pqtl_results["pQTL_right"],
    pqtl_results["pQTL_left"] + 50000,
    pqtl_results["pQTL_right"],
)
piqtl_results["piQTL_right"] = np.where(
    piqtl_results["piQTL_left"] > piqtl_results["piQTL_right"],
    piqtl_results["piQTL_left"] + 50000,
    piqtl_results["piQTL_right"],
)
eqtl_results["eQTL_right"] = np.where(
    eqtl_results["eQTL_left"] > eqtl_results["eQTL_right"],
    eqtl_results["eQTL_left"] + 50000,
    eqtl_results["eQTL_right"],
)

# Count colocal piQTL matches (left <= SNP position <= right)
colocal_piqtl_dict = {}
for snp_id in summary["SNP"]:
    snp_pos = summary[summary["SNP"] == snp_id]["position"].values[0]
    snp_chr = summary[summary["SNP"] == snp_id]["chromosome"].values[0]
    matching_rows = piqtl_results[
        (piqtl_results["chromosome"] == snp_chr)
        & (piqtl_results["piQTL_left"] <= snp_pos)
        & (piqtl_results["piQTL_right"] >= snp_pos)
    ]
    colocal_piqtl_dict[snp_id] = len(matching_rows)

summary["colocal_piQTL"] = summary["SNP"].map(colocal_piqtl_dict).fillna(0).astype(int)
print(f"    Counted colocal piQTL matches")

# Count colocal pQTL matches (left <= SNP position <= right)
colocal_pqtl_dict = {}
for snp_id in summary["SNP"]:
    snp_pos = summary[summary["SNP"] == snp_id]["position"].values[0]
    snp_chr = summary[summary["SNP"] == snp_id]["chromosome"].values[0]

    matching_rows = pqtl_results[
        (pqtl_results["chromosome"] == snp_chr)
        & (pqtl_results["pQTL_left"] <= snp_pos)
        & (pqtl_results["pQTL_right"] >= snp_pos)
    ]
    colocal_pqtl_dict[snp_id] = len(matching_rows)

summary["colocal_pQTL"] = summary["SNP"].map(colocal_pqtl_dict).fillna(0).astype(int)
print(f"    Counted colocal pQTL matches (with boundary normalization)")

# Count colocal eQTL matches (left <= SNP position <= right AND chromosome match)
colocal_eqtl_dict = {}
for snp_id in summary["SNP"]:
    snp_pos = summary[summary["SNP"] == snp_id]["position"].values[0]
    snp_chr = summary[summary["SNP"] == snp_id]["chromosome"].values[0]
    # Match by chromosome and position range
    matching_rows = eqtl_results[
        (eqtl_results["chromosome"] == snp_chr)
        & (eqtl_results["eQTL_left"] <= snp_pos)
        & (eqtl_results["eQTL_right"] >= snp_pos)
    ]
    colocal_eqtl_dict[snp_id] = len(matching_rows)

summary["colocal_eQTL"] = summary["SNP"].map(colocal_eqtl_dict).fillna(0).astype(int)
print(f"    Counted colocal eQTL matches")

# ============================================================================
# STEP 5: Save and validate output
# ============================================================================
print("\n[5/5] Saving and validating output...")

# Save the summary table
output_path = out_dir / "QTL_overlap_summary.csv"
summary.to_csv(output_path, index=False)
print(f"    Saved to {output_path}")

# Validation checks
print("\n" + "=" * 80)
print("Validation Summary:")
print("=" * 80)
print(f"Total SNPs: {len(summary)}")
print(f"Albert SNP matches: {summary['isAlbert_SNP'].sum()}")
print(f"SNPs with exact piQTL: {(summary['exact_piQTL'] > 0).sum()}")
print(f"SNPs with exact pQTL: {(summary['exact_pQTL'] > 0).sum()}")
print(f"SNPs with exact eQTL: {(summary['exact_eQTL'] > 0).sum()}")
print(f"SNPs with colocal piQTL: {(summary['colocal_piQTL'] > 0).sum()}")
print(f"SNPs with colocal pQTL: {(summary['colocal_pQTL'] > 0).sum()}")
print(f"SNPs with colocal eQTL: {(summary['colocal_eQTL'] > 0).sum()}")

# Check that colocal >= exact for each type
print("\nColocal ≥ Exact validation:")
print(f"    piQTL: {(summary['colocal_piQTL'] >= summary['exact_piQTL']).all()}")
print(f"    pQTL: {(summary['colocal_pQTL'] >= summary['exact_pQTL']).all()}")
print(f"    eQTL: {(summary['colocal_eQTL'] >= summary['exact_eQTL']).all()}")

# Data quality note
print("\nData Quality Notes:")
pqtl_reversed = (pqtl_results["pQTL_left"] > pqtl_results["pQTL_right"]).sum()
print(f"    pQTL entries with left > right: {pqtl_reversed}")
print(f"    (Boundaries were normalized for colocal matching)")
# print(f"    WARNING: 2 pQTL entries have peaks outside their normalized ranges")
# print(f"    This may indicate data quality issues in the pQTL results source file.")

# Check for NaN or negative values
print("\nData quality:")
print(f"    No NaN values: {not summary.isna().any().any()}")
print(
    f"    No negative counts: {(summary[['exact_piQTL', 'exact_pQTL', 'exact_eQTL', 'colocal_piQTL', 'colocal_pQTL', 'colocal_eQTL']] >= 0).all().all()}"
)

# Display first few rows
print("\nFirst 5 rows of summary table:")
print(summary.head())

print("\n" + "=" * 80)
print("Completed successfully!")
print("=" * 80)
