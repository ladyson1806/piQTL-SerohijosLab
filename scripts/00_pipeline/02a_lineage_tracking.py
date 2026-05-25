"""
02a_lineage_tracking.py
============================
Computes log2 ratio fitness scores from barcode read count data and generates
lineage tracking visualizations for protein-protein interaction (PPI) pairs
across drug and methotrexate (MTX) conditions.

Workflow
--------
For each PPI x DRUG combination:
  1. Read count data (strains x samples) are loaded once per PPI.
  2. Pseudocounts (+1) are added and counts are RPM-normalized.
  3. Log2(RPM_T96 / RPM_T0) fitness is computed for each strain, per
     replicate (A, B) and MTX condition (MTX, noMTX).
  4. Lineage tracking plots (RPM over time, log scale) and KDE fitness
     distributions are generated for each replicate x MTX condition.
  5. Replicate A vs B fitness scores are reference-normalized (subtracted
     mean of reference strains) and compared via Spearman correlation.
  6. Per-drug panels are assembled into a single vertical PNG per PPI x DRUG.
  7. A master fitness table and correlation statistics table are exported.

Inputs
------
../../data/pipeline/PPI_reference_barcodes.csv
    CSV with at least a 'PPI' column listing PPI identifiers to process.

../../results/01_barcode_count/per_PPI/{PPI}_read_number.csv
    Wide-format read count matrix for a given PPI.
    Rows: strains. Columns include:
      - strain_id
      - T0_noMTX_noDrug.{PPI}          : baseline read counts
      - {REP}_{MTX}_{DRUG}.{PPI}       : per-condition read counts
    where REP ∈ {A, B}, MTX ∈ {MTX, noMTX}, DRUG ∈ DRUGS constant.

Outputs
-------
../../results/02_lineage_tracking/{PPI}_{DRUG}.png
    Vertical panel per PPI x DRUG containing:
      - Lineage tracking + KDE plot for Replicate A (MTX and noMTX rows)
      - Lineage tracking + KDE plot for Replicate B (MTX and noMTX rows)
      - Replicate A vs B fitness correlation scatter (MTX and noMTX columns)

../../results/02_lineage_tracking/A_vs_B_fitness.csv
    Master table of reference-normalized fitness scores across all
    PPI x DRUG x MTX conditions, with columns:
      strain_id, logratio_Fitness_A, logratio_Fitness_B, MTX_condition

../../results/02_lineage_tracking/lineage_tracking_stats.csv
    Spearman correlation between replicate A and B fitness per condition,
    with columns: Condition, Correlation, p-value

Constants
---------
DRUGS       : list of str  — drug treatment labels
REPS        : list of str  — biological replicate labels ['A', 'B']
MTX_COND    : list of str  — methotrexate condition labels ['MTX', 'noMTX']
REF_STRAINS : list of str  — parent strain IDs used as fitness references
REF_COLORS  : list of str  — matplotlib colors assigned to reference strains in plots
REF_LABELS  : list of str  — display labels for reference strains in legends

Dependencies
------------
cv2, multiprocessing, os, signal, numpy, pandas,
seaborn, scipy, matplotlib, tqdm

Usage
-----
    python 02a_lineage_tracking.py
"""

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import cv2
import multiprocessing
import os
import signal
import numpy as np
import pandas as pd
import seaborn as sns
import scipy.stats as stats
import matplotlib
import matplotlib.pyplot as plt

from io import BytesIO
from tqdm import tqdm
from matplotlib.lines import Line2D

# ── Global plot style ────────────────────────────────────────────────────────
plt.rcParams.update({
    "figure.figsize":    [8, 6],
    "figure.autolayout": True,
    "font.family":       "sans-serif",
    "font.sans-serif":   "Helvetica",
    "font.size":         7,
})

DRUGS    = ['noDrug', '5-FC', 'Fluconazole', 'Metformin', 'Trifluoperazine']
REPS     = ['A', 'B']
MTX_COND = ['MTX', 'noMTX']
REF_STRAINS  = ['43', '599']
REF_COLORS   = ['paleturquoise', 'darkturquoise', 'darkcyan', 'yellow', 'orange', 'brown']
REF_LABELS   = ['43', '43_ref1', '43_ref2', '599', '599_ref1', '599_ref2']

# ── Helpers ──────────────────────────────────────────────────────────────────

def init_worker(tqdm_lock=None):
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    if tqdm_lock is not None:
        tqdm.set_lock(tqdm_lock)


def _read_img_from_buffer(buf: BytesIO) -> np.ndarray:
    """Decode a BytesIO PNG buffer into an OpenCV image array."""
    buf.seek(0)
    arr = np.frombuffer(buf.read(), dtype=np.uint8)
    img = cv2.imdecode(arr, cv2.IMREAD_COLOR)
    if img is None:
        raise ValueError("Failed to decode image from buffer.")
    return img


def vconcat_resize(img_list, interpolation=cv2.INTER_CUBIC):
    w_min = min(img.shape[1] for img in img_list)
    resized = [
        cv2.resize(img, (w_min, int(img.shape[0] * w_min / img.shape[1])), interpolation=interpolation)
        for img in img_list
    ]
    return cv2.vconcat(resized)


def hconcat_resize(img_list, interpolation=cv2.INTER_CUBIC):
    h_min = min(img.shape[0] for img in img_list)
    resized = [
        cv2.resize(img, (int(img.shape[1] * h_min / img.shape[0]), h_min), interpolation=interpolation)
        for img in img_list
    ]
    return cv2.hconcat(resized)


# ── Core computation ─────────────────────────────────────────────────────────

def compute_fitness(PPI_table: pd.DataFrame, PPI: str, REP: str, MTX: str, DRUG: str) -> pd.DataFrame:
    """
    Compute RPM-normalized log2 fitness for one REP x MTX x DRUG condition.
    Returns a tidy DataFrame with columns:
        strain_id, T0_RPM, T96_RPM, logratio_Fitness, Time (long format)
    and a wide DataFrame suitable for replicate merging.
    """
    t96_col = f'{REP}_{MTX}_{DRUG}.{PPI}'
    t0_col  = f'T0_noMTX_noDrug.{PPI}'

    R = PPI_table[['strain_id', t0_col, t96_col]].copy()
    R[t96_col] = R[t96_col] + 1
    R['T0_raw'] = R[t0_col] + 1

    rpm_t96_col = f'{t96_col}_RPM'
    rpm_t0_col  = f'T0.{PPI}_RPM'

    R[rpm_t96_col] = (R[t96_col] / R[t96_col].sum()) * 1e6
    R[rpm_t0_col]  = (R['T0_raw'] / R['T0_raw'].sum()) * 1e6
    R['logratio_Fitness'] = np.log2(R[rpm_t96_col] / R[rpm_t0_col])

    # Long format for lineage tracking plots
    t0_long  = R[['strain_id', rpm_t0_col,  'logratio_Fitness']].rename(columns={rpm_t0_col: 'RPM'})
    t0_long['Time'] = 0
    t96_long = R[['strain_id', rpm_t96_col, 'logratio_Fitness']].rename(columns={rpm_t96_col: 'RPM'})
    t96_long['Time'] = 96
    scatter_table = pd.concat([t0_long, t96_long], ignore_index=True)

    # Wide format for replicate correlation
    wide = R[['strain_id', rpm_t0_col, 'logratio_Fitness']].copy()

    return scatter_table, wide, rpm_t0_col


# ── Plotting helpers ─────────────────────────────────────────────────────────

def _make_legend_lines():
    return [Line2D([0], [0], color=c, linestyle='--', lw=2) for c in REF_COLORS]


def plot_lineage_rep(ax_lt, ax_kde, scatter_table: pd.DataFrame, title: str):
    """Plot lineage tracking lines + KDE fitness distribution for one REP x MTX."""
    segregants     = scatter_table[~scatter_table['strain_id'].str.contains('ref', case=False)]
    ref_segregants = scatter_table[
        scatter_table['strain_id'].str.contains('ref', case=False) |
        scatter_table['strain_id'].isin(REF_STRAINS)
    ].reset_index(drop=True)

    sns.lineplot(data=segregants, x='Time', y='RPM',
                 hue='logratio_Fitness', legend=False, ax=ax_lt)
    sns.lineplot(data=ref_segregants, x='Time', y='RPM',
                 hue='strain_id', palette=REF_COLORS,
                 linestyle='dashed', legend=False, ax=ax_lt)
    ax_lt.set_yscale('log')
    ax_lt.set_title(title)

    sns.kdeplot(data=scatter_table['logratio_Fitness'], ax=ax_kde)
    ax_kde.set_title(title)

    for k, ref in enumerate(np.unique(ref_segregants['strain_id'])):
        ref_fit = ref_segregants.loc[ref_segregants['strain_id'] == ref, 'logratio_Fitness'].iloc[0]
        ax_kde.axvline(x=ref_fit, linestyle='dashed', color=REF_COLORS[k])

    return ref_segregants  # returned so caller can build legend


# ── Main per-PPI/DRUG function ───────────────────────────────────────────────

def fitness_comparison(PPI: str, DRUG: str, res: dict):
    """
    For a given PPI x DRUG:
      - Compute fitness for each REP x MTX combination
      - Plot lineage tracking + KDE panels (MTX / noMTX) per replicate
      - Plot replicate A vs B fitness correlation
      - Assemble vertical panel and write to results folder
    Returns (fitness_table, res, output_path).
    """
    # Load CSV once
    PPI_table = pd.read_csv(f'../../results/01_barcode_count/per_PPI/{PPI}_read_number.csv')

    # Precompute fitness for all 4 conditions: REP x MTX
    fitness_data = {}  # key: (REP, MTX)
    for REP in REPS:
        for MTX in MTX_COND:
            scatter_table, wide, t0_rpm_col = compute_fitness(PPI_table, PPI, REP, MTX, DRUG)
            fitness_data[(REP, MTX)] = {
                'scatter': scatter_table,
                'wide':    wide,
                't0_col':  t0_rpm_col,
            }

    # ── Lineage tracking panels (one per REP, both MTX conditions side by side) ──
    lt_buffers = {}
    for REP in REPS:
        g, g_axes = plt.subplots(ncols=2, nrows=2, figsize=(12, 12))
        legend_lines = _make_legend_lines()
        g_axes[0, 1].legend(legend_lines, REF_LABELS, loc='center right', fontsize=8)
        g_axes[1, 1].legend(legend_lines, REF_LABELS, loc='center right', fontsize=8)

        for row, MTX in enumerate(MTX_COND):
            title = f'{PPI}_{DRUG}_{MTX} - Replicate {REP}'
            plot_lineage_rep(
                g_axes[row, 0], g_axes[row, 1],
                fitness_data[(REP, MTX)]['scatter'],
                title
            )

        buf = BytesIO()
        g.savefig(buf, format='png', dpi=300)
        plt.close(g)
        lt_buffers[REP] = buf

    # ── Replicate A vs B fitness correlation (MTX and noMTX separately) ──
    f, axes = plt.subplots(ncols=2, sharex=True, sharey=True, figsize=(12, 8))

    all_merged = []
    for i, MTX in enumerate(MTX_COND):
        wide_A = fitness_data[('A', MTX)]['wide'].rename(columns={'logratio_Fitness': 'logratio_Fitness_A'})
        wide_B = fitness_data[('B', MTX)]['wide'].rename(columns={'logratio_Fitness': 'logratio_Fitness_B'})
        t0_col = fitness_data[('A', MTX)]['t0_col']

        REP_MERGED = pd.merge(wide_A, wide_B, on=['strain_id', t0_col])
        REP_MERGED_REF = REP_MERGED[REP_MERGED['strain_id'].str.contains('_ref')]

        REP_MERGED['logratio_Fitness_A'] -= REP_MERGED_REF['logratio_Fitness_A'].mean()
        REP_MERGED['logratio_Fitness_B'] -= REP_MERGED_REF['logratio_Fitness_B'].mean()
        REP_MERGED['MTX_condition'] = f'{PPI}_{DRUG}_{MTX}'

        sns.scatterplot(data=REP_MERGED, x='logratio_Fitness_A', y='logratio_Fitness_B',
                        color='black', ax=axes[i])
        sns.scatterplot(data=REP_MERGED[REP_MERGED['strain_id'].str.contains('_ref')],
                        x='logratio_Fitness_A', y='logratio_Fitness_B',
                        color='red', ax=axes[i])

        corr, pval = stats.spearmanr(REP_MERGED['logratio_Fitness_A'], REP_MERGED['logratio_Fitness_B'])
        axes[i].set_xlabel('Fitness (Replicate A)')
        axes[i].set_ylabel('Fitness (Replicate B)')
        axes[i].annotate(
            f"Spearman r: {corr:.3f}\np-val: {pval:.2e}",
            xy=[0, -60], xycoords="axes points"
        )
        axes[i].set_title(f'{PPI} under {DRUG} ({MTX})'.replace('_', ':'))

        res[f'{PPI}_{DRUG}_{MTX}'] = {'corr': corr, 'pval': pval}
        all_merged.append(REP_MERGED[['strain_id', 'logratio_Fitness_A', 'logratio_Fitness_B', 'MTX_condition']])

    corr_buf = BytesIO()
    f.savefig(corr_buf, format='png', dpi=300)
    plt.close(f)

    # ── Assemble vertical panel ──────────────────────────────────────────────
    imgs = [
        _read_img_from_buffer(lt_buffers['A']),
        _read_img_from_buffer(lt_buffers['B']),
        _read_img_from_buffer(corr_buf),
    ]
    panel = vconcat_resize(imgs)
    out_path = f'../../results/02_lineage_tracking/{PPI}_{DRUG}.png'
    cv2.imwrite(out_path, panel)

    return pd.concat(all_merged, ignore_index=True), res, out_path


# ── Entry point ───────────────────────────────────────────────

# Load PPI list and set up results folder
PPI_list = pd.read_csv('../../data/pipeline/PPI_reference_barcodes.csv')['PPI']
results_folder = '../../results'

# Create results subfolder if it doesn't exist
for fld in [
    os.path.join(results_folder, '02_lineage_tracking')
]:
    os.makedirs(fld, exist_ok=True)

# For testing, we can limit to a subset of PPIs. In production, use the full list from the CSV.
# PPI_list = ['ALO1_ADE17', 'ERG11_MID2']
ALL_TABLE = []
res = {}
merge_args = []

# Process each PPI x DRUG combination and collect fitness tables and stats
for PPI in tqdm(PPI_list, desc="PPI"):
    ALL_DRUGS = []
    print(f'\nLineage Tracking for {PPI} in progress')
    for DRUG in tqdm(DRUGS, desc=f"  {PPI} drugs", leave=False):
        TABLE, res, drug_image = fitness_comparison(PPI, DRUG, res)
        ALL_DRUGS.append(drug_image)
        ALL_TABLE.append(TABLE)

# Save master fitness table
pd.concat(ALL_TABLE, ignore_index=True).to_csv('../../results/02_lineage_tracking/A_vs_B_fitness.csv', index=False)

# Save correlation stats
stats_rows = [
    {'Condition': k, 'Correlation': v['corr'], 'p-value': v['pval']}
    for k, v in res.items()
]
pd.DataFrame(stats_rows).to_csv('../../results/02_lineage_tracking/lineage_tracking_stats.csv', index=False)