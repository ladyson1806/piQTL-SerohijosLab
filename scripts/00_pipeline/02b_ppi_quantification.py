# =============================================================================
# PPI Fitness Estimation Pipeline
# =============================================================================
# Steps:
#   1. Split barcode counts per PPI
#   2. Compute log-ratio fitness per PPI x Drug x Replicate
#   3. Merge replicates and drugs per PPI
#   4. Normalize by subtracting reference strain mean (v1)
#   5. Normalize by subtracting reference strain mean + noPPI baseline (v2)
# =============================================================================

import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)

import os
import numpy as np
import pandas as pd
from tqdm import tqdm
from functools import reduce

# =============================================================================
# Constants
# =============================================================================

DRUGS = ["noDrug", "5-FC", "Fluconazole", "Metformin", "Trifluoperazine"]
MTX_CONDITIONS = ["noMTX", "MTX"]
REPLICATES = ["A", "B"]

LOGRATIO_CUSTOM_COLS = [
    "strain_id",
    "noMTX_noDrug_avg_logratio_Fitness",   "MTX_noDrug_avg_logratio_Fitness",
    "noMTX_5-FC_avg_logratio_Fitness",     "MTX_5-FC_avg_logratio_Fitness",
    "noMTX_Fluconazole_avg_logratio_Fitness", "MTX_Fluconazole_avg_logratio_Fitness",
    "noMTX_Metformin_avg_logratio_Fitness",   "MTX_Metformin_avg_logratio_Fitness",
    "noMTX_Trifluoperazine_avg_logratio_Fitness", "MTX_Trifluoperazine_avg_logratio_Fitness",
]

# =============================================================================
# Helper functions
# =============================================================================

def merge_on_strain(tables):
    """Outer-merge a list of DataFrames on 'strain_id'."""
    return reduce(lambda l, r: pd.merge(l, r, on="strain_id", how="outer"), tables)


def build_metadata_table(columns):
    """
    Parse column names of the form '{PPI1}_{PPI2}_{MTX}_{DRUG}_...'
    into a metadata DataFrame with columns [label, PPI, Drug, MTX].
    """
    records = []
    for col in columns:
        if col == "TAXA":
            continue
        parts = col.split("_")
        records.append({
            "label": col,
            "PPI":   f"{parts[0]}_{parts[1]}",
            "MTX":   parts[2],
            "Drug":  parts[3],
        })
    return pd.DataFrame(records, columns=["label", "PPI", "Drug", "MTX"])


def load_logratio_tables(PPI_list, results_folder):
    """Load per-PPI logratio fitness CSVs with fixed column names."""
    return [
        pd.read_csv(
            f"{results_folder}/03_ppi_estimation/logratio/before_downsampling/{PPI}_logratio_fitness.csv",
            header=None, names=LOGRATIO_CUSTOM_COLS,
        ).drop(0).reset_index(drop=True)
        for PPI in PPI_list
    ]


def cast_fitness_cols_to_float(df):
    """Cast all columns except 'strain_id' to float, in place."""
    for col in df.columns:
        if col != "strain_id":
            df[col] = df[col].astype(float)


def filter_ref_strains(df, id_col="strain_id"):
    """Remove reference strains and rename the ID column to 'TAXA'."""
    return (
        df[~df[id_col].str.contains("_ref")]
        .rename(columns={id_col: "TAXA"})
    )


def save_for_eqtl_and_rmvp(df, genotype, results_folder, tag):
    """Export a normalized PPI table for eQTL matrix and rMVP formats."""
    base = f"{results_folder}/03_ppi_estimation/logratio"
    pivot = df.rename(columns={"TAXA": "Condition"}).set_index("Condition").T
    strain_cols = genotype.columns[1:]

    pivot[strain_cols].to_csv(
        f"{base}/all_PPI_logratio_fitness_{tag}_for_eQTL_matrix.csv"
    )
    (
        pivot[strain_cols].T
        .reset_index()
        .rename(columns={"Condition": "TAXA"})
        .to_csv(f"{base}/all_PPI_logratio_fitness_{tag}_for_rMVP.csv", index=False)
    )

# =============================================================================
# Step 1 — Split barcode counts per PPI
# =============================================================================

def reads_per_ppi(barcode_count_folder, PPI_list, results_folder):
    """Merge per-condition barcode count CSVs and split by PPI."""
    tables = [
        pd.read_csv(os.path.join(barcode_count_folder, f))
        for f in os.listdir(barcode_count_folder) if f.endswith(".csv")
    ]
    ALL = merge_on_strain(tables).replace(np.nan, 0)

    for PPI in sorted(PPI_list):
        cols = ["strain_id"] + sorted(col for col in ALL.columns if PPI in col)
        ALL[cols].to_csv(
            f"{results_folder}/01_barcode_count/per_PPI/{PPI}_read_number.csv",
            index=False,
        )

# =============================================================================
# Step 2 — Compute log-ratio fitness for one PPI x Drug
# =============================================================================

def logratio_fitness(PPI, DRUG):
    """
    For a given PPI and drug, compute RPM-normalized log-ratio fitness
    for each MTX condition, averaged across replicates A and B.
    Returns a DataFrame with columns [strain_id, {PPI}_noMTX_{DRUG}_avg_logratio_Fitness,
                                                  {PPI}_MTX_{DRUG}_avg_logratio_Fitness].
    """
    PPI_table = pd.read_csv(
        f"../../results/01_barcode_count/per_PPI/{PPI}_read_number.csv"
    )

    all_mtx = []
    for MTX in MTX_CONDITIONS:
        all_rep = []
        for REP in REPLICATES:
            R = PPI_table[["strain_id", f"T0_noMTX_noDrug.{PPI}", f"{REP}_{MTX}_{DRUG}.{PPI}"]].copy()

            # Pseudocount
            R[f"{REP}_{MTX}_{DRUG}.{PPI}"] += 1
            R[f"T0.{PPI}"] = R[f"T0_noMTX_noDrug.{PPI}"] + 1

            # RPM normalization
            R[f"{REP}_{MTX}_{DRUG}.{PPI}_RPM"] = (
                R[f"{REP}_{MTX}_{DRUG}.{PPI}"] / R[f"{REP}_{MTX}_{DRUG}.{PPI}"].sum()
            ) * 1e6
            R[f"T0.{PPI}_RPM"] = (R[f"T0.{PPI}"] / R[f"T0.{PPI}"].sum()) * 1e6

            # Log-ratio fitness
            fitness_col = f"{PPI}_{MTX}_{DRUG}_logratio_Fitness"
            R[fitness_col] = np.log2(
                R[f"{REP}_{MTX}_{DRUG}.{PPI}_RPM"] / R[f"T0.{PPI}_RPM"]
            )

            keep_cols = ["strain_id"] + [c for c in R.columns if "RPM" in c] + [fitness_col]
            all_rep.append(R[keep_cols])

        # Average fitness across replicates
        merged = pd.merge(
            all_rep[0], all_rep[1],
            on=["strain_id", f"T0.{PPI}_RPM"],
            suffixes=("_A", "_B"),
        )
        merged[f"{PPI}_{MTX}_{DRUG}_avg_logratio_Fitness"] = (
            merged[f"{PPI}_{MTX}_{DRUG}_logratio_Fitness_A"]
            + merged[f"{PPI}_{MTX}_{DRUG}_logratio_Fitness_B"]
        ) / 2
        all_mtx.append(merged[["strain_id", f"{PPI}_{MTX}_{DRUG}_avg_logratio_Fitness"]])

    return pd.merge(all_mtx[0], all_mtx[1], on="strain_id")

# =============================================================================
# Step 3 — Run logratio fitness across all PPIs and drugs
# =============================================================================

def run_logratio_fitness_per_ppi(PPI_list, results_folder):
    """Compute and save logratio fitness for every PPI across all drugs."""
    for PPI in tqdm(PPI_list, desc="Computing logratio fitness"):
        drug_tables = [logratio_fitness(PPI, DRUG) for DRUG in DRUGS]
        FINAL = merge_on_strain(drug_tables)

        logratio_cols = ["strain_id"] + [c for c in FINAL.columns if "_logratio_Fitness" in c]
        output = os.path.join(
            results_folder,
            f"03_ppi_estimation/logratio/before_downsampling/{PPI}_logratio_fitness.csv",
        )
        FINAL[logratio_cols].to_csv(output, index=False)

# =============================================================================
# Step 4 — Merge all PPIs before normalization
# =============================================================================

def ppi_quantification_before_normalization(PPI_list, results_folder):
    """Concatenate per-PPI logratio tables, remove ref strains, save + metadata."""
    ppi_tables = [
        pd.read_csv(
            f"{results_folder}/03_ppi_estimation/logratio/before_downsampling/{PPI}_logratio_fitness.csv"
        )
        for PPI in PPI_list
    ]
    all_logratio = filter_ref_strains(merge_on_strain(ppi_tables))

    all_logratio.to_csv(
        os.path.join(results_folder, "03_ppi_estimation/logratio/all_PPI_logratio_fitness.csv"),
        index=False,
    )

    metadata = build_metadata_table(all_logratio.columns)
    metadata.to_csv(
        os.path.join(results_folder, "03_ppi_estimation/logratio/all_PPI_logratio_metadata.csv"),
        index=False,
    )

# =============================================================================
# Step 5a — Normalize by subtracting reference strain mean
# =============================================================================

def ppi_quantification_after_normalization(PPI_list, results_folder, genotype):
    """Subtract per-PPI reference strain mean; export eQTL and rMVP matrices."""
    ppi_tables = load_logratio_tables(PPI_list, results_folder)

    all_ppi = []
    for i, PPI in enumerate(PPI_list):
        df = ppi_tables[i]
        cast_fitness_cols_to_float(df)

        ref_mean = df[df["strain_id"].str.contains("_ref")].mean(numeric_only=True)
        fitness_cols = [c for c in df.columns if c != "strain_id"]

        delta = df.copy()
        for col in fitness_cols:
            delta[f"{col}_minus_ref"] = delta[col] - ref_mean[col]

        minus_ref_cols = ["strain_id"] + [c for c in delta.columns if c.endswith("_minus_ref")]
        delta = delta[minus_ref_cols]
        delta.columns = ["strain_id"] + [f"{PPI}_{c}" for c in delta.columns[1:]]
        all_ppi.append(delta)

    normalized = filter_ref_strains(merge_on_strain(all_ppi))
    save_for_eqtl_and_rmvp(
        normalized, genotype, results_folder,
        tag="before_downsampling_minus_ref",
    )

    metadata = build_metadata_table(normalized.columns)
    metadata.to_csv(
        f"{results_folder}/03_ppi_estimation/logratio/all_PPI_logratio_minus_ref_metadata.csv",
        index=False,
    )

# =============================================================================
# Step 5b — Normalize by subtracting ref mean AND noPPI baseline
# =============================================================================

def ppi_quantification_after_normalization_v2(PPI_list, results_folder, genotype):
    """Subtract ref strain mean then subtract the noPPI reference table."""
    ppi_tables = load_logratio_tables(PPI_list, results_folder)

    # Last entry in PPI_list is the noPPI control
    noPPI_table = ppi_tables[-1].copy()
    cast_fitness_cols_to_float(noPPI_table)

    all_ppi = []
    for i, PPI in enumerate(PPI_list[:-1]):  # exclude noPPI control
        df = ppi_tables[i].copy()
        cast_fitness_cols_to_float(df)

        ref_mean = df[df["strain_id"].str.contains("_ref")].mean(numeric_only=True)
        for col in ref_mean.index:
            df[col] = df[col] - ref_mean[col]

        delta = (
            df.set_index("strain_id")
            .subtract(noPPI_table.set_index("strain_id"))
            .reset_index()
        )
        delta.columns = ["strain_id"] + [f"{PPI}_{c}" for c in delta.columns[1:]]
        all_ppi.append(delta)

    normalized = filter_ref_strains(merge_on_strain(all_ppi))
    save_for_eqtl_and_rmvp(
        normalized, genotype, results_folder,
        tag="before_downsampling_minus_ref_noPPI-labelled",
    )

    metadata = build_metadata_table(normalized.columns)
    metadata.to_csv(
        f"{results_folder}/03_ppi_estimation/logratio/all_PPI_logratio_minus_ref_noPPI-labelled_metadata.csv",
        index=False,
    )

# =============================================================================
# Entry point
# =============================================================================

if __name__ == "__main__":
    genotype   = pd.read_csv("../../data/genotype_information/piQTL_genotype_matrix_dec2022.txt")
    PPI_list   = pd.read_csv("../../data/pipeline/PPI_reference_barcodes.csv")["PPI"].tolist()

    results_folder       = "../../results/"
    per_condition_folder = os.path.join(results_folder, "01_barcode_count/per_condition")

    # Create output directories
    for folder in [
        "03_ppi_estimation",
        "03_ppi_estimation/logratio",
        "03_ppi_estimation/logratio/before_downsampling",
    ]:
        os.makedirs(os.path.join(results_folder, folder), exist_ok=True)

    print("Step 1 — Splitting barcode counts per PPI ...")
    reads_per_ppi(per_condition_folder, PPI_list, results_folder)

    print("Step 2-3 — Estimating fitness from barcode frequency ...")
    run_logratio_fitness_per_ppi(PPI_list, results_folder)

    print("Step 4 — Merging replicates across all PPIs ...")
    ppi_quantification_before_normalization(PPI_list, results_folder)

    print("Step 5a — Normalizing by reference strain mean ...")
    ppi_quantification_after_normalization(PPI_list, results_folder, genotype)

    print("Step 5b — Normalizing by reference strain mean + noPPI baseline ...")
    ppi_quantification_after_normalization_v2(PPI_list, results_folder, genotype)