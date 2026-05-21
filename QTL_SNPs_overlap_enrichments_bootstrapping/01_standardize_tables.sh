#!/bin/sh

cd $(dirname ${0})

mkdir -p out/standardized_tables


# Standardize the piQTL SNP annotation table
STANDARDIZED_PIQTL_SNP_ANNOTATION="out/standardized_tables/piQTL_SNP_annotation.csv"
if [ ! -e ${STANDARDIZED_PIQTL_SNP_ANNOTATION} ]; then
    python src/standardize_snp_table.py \
        --input_file data/piQTL_SNP_annotation.csv \
        --output_file ${STANDARDIZED_PIQTL_SNP_ANNOTATION} \
        --blacklist data/masked_SNPs.txt
fi


# Standardize the pQTL results table
STANDARDIZED_PQTL_RESULTS="out/standardized_tables/pQTL_results.csv"
if [ ! -e ${STANDARDIZED_PQTL_RESULTS} ]; then
    python src/standardize_pqtl_results.py \
        --input_file data/pQTL_results.csv \
        --output_file ${STANDARDIZED_PQTL_RESULTS}
fi


# Standardize the eQTL results table
STANDARDIZED_EQTL_RESULTS="out/standardized_tables/eQTL_results.csv"
if [ ! -e ${STANDARDIZED_EQTL_RESULTS} ]; then
    python src/standardize_eqtl_results.py \
        --input_file data/eQTL_results.csv \
        --output_file ${STANDARDIZED_EQTL_RESULTS}
fi


# Standardize the SNP annotation table with whitelist (actually detected piQTL SNPS).
STANDARDIZED_PIQTL_SNP_ANNOTATION_WITH_WHITELIST="out/standardized_tables/piQTL_SNP_annotation_with_whitelist.csv"
if [ ! -e ${STANDARDIZED_PIQTL_SNP_ANNOTATION_WITH_WHITELIST} ]; then
    python src/standardize_snp_table.py \
        --input_file data/piQTL_SNP_annotation.csv \
        --output_file ${STANDARDIZED_PIQTL_SNP_ANNOTATION_WITH_WHITELIST} \
        --whitelist data/actual_piQTL_SNPs.txt
fi


# Standardize the SNP annotation table without any filtering (for the purpose of checking the effect of filtering).
STANDARDIZED_PIQTL_SNP_ANNOTATION_NO_FILTERING="out/standardized_tables/piQTL_SNP_annotation_no_filtering.csv"
if [ ! -e ${STANDARDIZED_PIQTL_SNP_ANNOTATION_NO_FILTERING} ]; then
    python src/standardize_snp_table.py \
        --input_file data/piQTL_SNP_annotation.csv \
        --output_file ${STANDARDIZED_PIQTL_SNP_ANNOTATION_NO_FILTERING}
fi


# Standardize the SNP annotation table for the eQTLs for confirming the shared SNPs between piQTLs and eQTLs.
STANDARDIZED_EQTL_SNP_ANNOTATION="out/standardized_tables/eQTL_SNP_annotation.csv"
if [ ! -e ${STANDARDIZED_EQTL_SNP_ANNOTATION} ]; then
    python src/standardize_eqtl_snp_table.py \
        --input_file data/eQTL_SNP_annotation.csv \
        --output_file ${STANDARDIZED_EQTL_SNP_ANNOTATION}
fi


# Extract the shared SNPs between piQTLs and eQTLs.
STANDARDIZED_SHARED_SNP_ANNOTATION="out/standardized_tables/piQTL_SNP_annotation_shared_with_eQTLs.csv"
if [ ! -e ${STANDARDIZED_SHARED_SNP_ANNOTATION} ]; then
    python src/extract_shared_snps.py \
        --piqtl_snp_annotation ${STANDARDIZED_PIQTL_SNP_ANNOTATION} \
        --eqtl_snp_annotation ${STANDARDIZED_EQTL_SNP_ANNOTATION} \
        --output_file ${STANDARDIZED_SHARED_SNP_ANNOTATION}
fi
