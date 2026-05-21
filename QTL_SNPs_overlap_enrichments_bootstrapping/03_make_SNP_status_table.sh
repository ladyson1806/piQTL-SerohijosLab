#!/bin/sh

cd $(dirname ${0})

mkdir -p out/overlap_status

# Input files
PIQTL_TABLE="out/standardized_tables/piQTL_SNP_annotation.csv"
PIQTL_TABLE_NO_FILTERING="out/standardized_tables/piQTL_SNP_annotation_no_filtering.csv"
PQTL_TABLE="out/standardized_tables/pQTL_results.csv"
EQTL_TABLE="out/standardized_tables/eQTL_results.csv"
ACTUAL_SNP_FILE="data/actual_piQTL_SNPs.txt"

OVERLAP_STATUS_TABLE="out/overlap_status/piQTL_overlap_status.csv"
if [ ! -f ${OVERLAP_STATUS_TABLE} ]; then
    echo "Generate overlap status table...."
    python src/build_piQTL_overlap_status.py \
        --piqtl ${PIQTL_TABLE} \
        --pqtl ${PQTL_TABLE} \
        --eqtl ${EQTL_TABLE} \
        --actual ${ACTUAL_SNP_FILE} \
        --out ${OVERLAP_STATUS_TABLE}
fi


# Also generate the overlap status table without any filtering of piQTL SNPs, for the purpose of checking the effect of filtering.
OVERLAP_STATUS_TABLE_NO_FILTERING="out/overlap_status/piQTL_overlap_status_no_filtering.csv"
if [ ! -f ${OVERLAP_STATUS_TABLE_NO_FILTERING} ]; then
    echo "Generate overlap status table without filtering...."
    python src/build_piQTL_overlap_status.py \
        --piqtl ${PIQTL_TABLE_NO_FILTERING} \
        --pqtl ${PQTL_TABLE} \
        --eqtl ${EQTL_TABLE} \
        --actual ${ACTUAL_SNP_FILE} \
        --out ${OVERLAP_STATUS_TABLE_NO_FILTERING}
fi
