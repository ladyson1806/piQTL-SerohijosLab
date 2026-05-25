#!/bin/sh

cd $(dirname ${0})

mkdir -p out/randomized_overlap_counts_exact
mkdir -p out/randomized_overlap_counts_colocal

# Randomized overlap analysis settings
N_ITERATIONS=10000
RANDOM_SEED=42
NUM_SNP=354

PIQTL_TABLE="out/standardized_tables/piQTL_SNP_annotation.csv"
PQTL_TABLE="out/standardized_tables/pQTL_results.csv"
EQTL_TABLE="out/standardized_tables/eQTL_results.csv"

python src/overlap_counts_exact.py \
	--n ${N_ITERATIONS} \
	--seed ${RANDOM_SEED} \
	--num_snp ${NUM_SNP} \
	--piqtl ${PIQTL_TABLE} \
	--pqtl ${PQTL_TABLE} \
	--eqtl ${EQTL_TABLE} \
	--outdir out/randomized_overlap_counts_exact

python src/overlap_counts_colocal.py \
	--n ${N_ITERATIONS} \
	--seed ${RANDOM_SEED} \
	--num_snp ${NUM_SNP} \
	--piqtl ${PIQTL_TABLE} \
	--pqtl ${PQTL_TABLE} \
	--eqtl ${EQTL_TABLE} \
	--outdir out/randomized_overlap_counts_colocal
