# piQTL — Serohijos Lab

**piQTL** (protein-protein interaction QTL) is a genome-wide mapping framework that links natural genetic variation in *Saccharomyces cerevisiae* to quantitative changes in the formation of protein-protein interaction (PPI) complexes, measured across multiple drug environments.

This repository contains the full computational pipeline, downstream analysis modules, and interactive visualization apps associated with the published manuscript.

> **Citation:** _[add DOI / full citation here]_
>
> Interactive data portal: [ladyson1806.github.io/piQTL-SerohijosLab](https://ladyson1806.github.io/piQTL-SerohijosLab/)

Zenodo: [![DOI](https://zenodo.org/badge/1238377380.svg)](https://doi.org/10.5281/zenodo.21069963)


---

## Overview

The piQTL approach uses a barcoded yeast strain collection and DHFR protein-fragment complementation assay (PCA) to quantify PPI fitness across a population of inbred lines. Sequencing barcode counts over time (lineage tracking) produces per-strain, per-condition fitness scores that serve as quantitative traits. Genome-wide association with SNP genotypes identifies **piQTLs** — loci whose allelic state affects PPI complex abundance.

Drug conditions screened: **noDrug**, **5-FC**, **Fluconazole**, **Metformin**, **Trifluoperazine** (each ± methotrexate / MTX).

---

## Repository Structure

```
scripts/
├── 00_pipeline/             # End-to-end pipeline from raw sequencing to piQTLs
├── 01_experimental_checks/  # Quality control and LD characterization
└── 02_downstream_analyses/  # Post-mapping comparative and enrichment analyses

shiny_apps/
├── piQTL_manhattan/         # Interactive Manhattan & QQ plots (R Shiny)
└── piQTL_genome_browser/    # Interactive genome browser for piQTL tracks (R Shiny)

data/                        # Reference annotations and genotype information
results/                     # Pipeline outputs
docs/                        # GitHub Pages source (Jekyll)
figures/                     # Manuscript figures
```

---

## Pipeline (`scripts/00_pipeline/`)

The main pipeline runs in numbered order:

| Step | Script | Description |
|------|--------|-------------|
| 00a | `00a_generate_genotype_matrix.py` | Build SNP genotype matrix from She & Jarosz 2018 phased genotypes |
| 00b | `00b_get_yeast_genome_annotations.ipynb` | Download and format SGD genome annotations |
| 00c | `00c_barcode_reference_library.py` | Assemble the barcode-to-strain reference library |
| 01a | `01a_fastq_concatenation.ipynb` | Concatenate raw FASTQ files per sample |
| 01b | `01b_fastq_preprocessing.py` | Trim and quality-filter reads |
| 01c | `01c_barcode_extraction.py` | Extract barcodes from preprocessed reads |
| 01d | `01d_barcode_mapping.py` | Map barcodes to the reference library |
| 02a | `02a_lineage_tracking.py` | Compute log₂-ratio fitness scores and lineage tracking plots |
| 02b | `02b_ppi_quantification.py` | Estimate per-PPI fitness across drug × MTX conditions |
| 03  | `03_run_eQTL_matrix.R` | Run genome-wide QTL scan with MatrixEQTL |
| 04  | `04_piQTLs_collection_and_analyses.ipynb` | Collect significant piQTLs and apply masking / LD clumping |

---

## Experimental Checks (`scripts/01_experimental_checks/`)

- **LD characterization** — `generate_LD_table.R`, `LD_clumping.ipynb`: compute pairwise LD across the genotype panel and define LD blocks.
- **Genotype heatmap** — `build_genotype_heatmap.py`: visualize strain-level genotype similarity.
- **Replicability** — `experimental_design_replicability.py`: assess concordance between biological replicates.

---

## Downstream Analyses (`scripts/02_downstream_analyses/`)

| Module | Description |
|--------|-------------|
| `comparing_piQTL_pQTL_eQTL/` | Exact-match and LD-colocalized overlap tables across piQTL, pQTL, and eQTL datasets |
| `QTL_mapping_on_genome/` | Genome-wide summary tables and visualization of QTL overlap counts |
| `QTL_SNPs_overlap_enrichments_bootstrapping/` | Randomized SNP sampling and empirical p-value calculation for overlap enrichments |
| `cis_trans_QTL_effects/` | Cis/trans effect-size processing from multiple source datasets (piQTL, natural isolates, human studies) |
| `cis_trans_QTL_effects_bootstrapping/` | Bootstrap simulations and summary-statistics comparison for cis/trans metrics |
| `PCAs/` | Principal component analyses of PPI fitness landscapes |
| `piQTL_plots/` | Manuscript-ready figure generation |
| `network_analyses/` | PPI similarity network construction and community analysis |
| `rnaseq_analysis/` | Differential expression analysis (DESeq2) for RNA-seq validation |

---

## Interactive Visualization Apps

Both Shiny apps are publicly hosted:

- **Manhattan & QQ plots** — [serohijos-piqtl.shinyapps.io/piQTL_manhattan](https://serohijos-piqtl.shinyapps.io/piQTL_manhattan/)  
  Select a PPI and drug condition to display genome-wide association results.

- **Genome Browser** — [serohijos-piqtl.shinyapps.io/piQTL_genome_browser](https://serohijos-piqtl.shinyapps.io/piQTL_genome_browser/)  
  Load piQTL tracks alongside LD blocks and non-coding RNA annotations. Search by gene name or genomic coordinates.

Source code for both apps lives in `shiny_apps/`.

---

## Environment

**OS:** Linux / macOS recommended.

**Python 3.x** — core pipeline and downstream analyses:
```
pandas, numpy, scipy, matplotlib, seaborn, plotly, tqdm, roman
```

**R** — QTL mapping and Shiny apps:
```
MatrixEQTL, DESeq2, stringr, readr, dplyr, tibble, shiny
```

Install Python dependencies with:
```bash
pip install pandas numpy scipy matplotlib seaborn plotly tqdm roman
```

---

## Data

- **Raw sequencing data:** available on NCBI BioProject _(link to be added)_.
- **Genotype matrix:** `data/genotype_information/piQTL_genotype_matrix_dec2022.txt` — derived from She & Jarosz 2018 (*S. cerevisiae* inbred panel).
- **Genome annotations:** `data/genome_annotations/` — SGD R64-3-1 ORF annotations, LD block definitions, and gene-by-strain tables.


---

## Notes

- This repository is a reproducibility resource. It does not provide biological interpretation of results; please refer to the manuscript for that.
- Each subdirectory under `scripts/02_downstream_analyses/` contains a local README describing inputs, outputs, and usage.
- Input data for downstream modules are expected in each module's local `data/` directory.
