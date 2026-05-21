#!/bin/sh

cd $(dirname ${0})

mkdir -p out/test_randomized_overlap_counts_exact
mkdir -p out/test_randomized_overlap_counts_colocal

# Randomized overlap analysis settings
N_ITERATIONS=10
RANDOM_SEED=1234
NUM_SNP=354

PIQTL_TABLE="out/standardized_tables/piQTL_SNP_annotation_with_whitelist.csv"
PQTL_TABLE="out/standardized_tables/pQTL_results.csv"
EQTL_TABLE="out/standardized_tables/eQTL_results.csv"

python src/overlap_counts_exact.py \
	--n ${N_ITERATIONS} \
	--seed ${RANDOM_SEED} \
	--num_snp ${NUM_SNP} \
	--piqtl ${PIQTL_TABLE} \
	--pqtl ${PQTL_TABLE} \
	--eqtl ${EQTL_TABLE} \
	--outdir out/test_randomized_overlap_counts_exact

python src/overlap_counts_colocal.py \
	--n ${N_ITERATIONS} \
	--seed ${RANDOM_SEED} \
	--num_snp ${NUM_SNP} \
	--piqtl ${PIQTL_TABLE} \
	--pqtl ${PQTL_TABLE} \
	--eqtl ${EQTL_TABLE} \
	--outdir out/test_randomized_overlap_counts_colocal
