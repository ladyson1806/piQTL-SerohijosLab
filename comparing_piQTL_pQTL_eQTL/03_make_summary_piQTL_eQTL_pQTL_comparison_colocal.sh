#!/bin/sh

cd $(dirname ${0})

mkdir -p out/summary_colocal


# Set template SNP table for storing pi/e/pQTLs.
MASTER_TABLE="out/summary/summary_master_SNP_tabl.csv"


# Summarize colocaled pQTLs and eQTLs grouped by piQTLs into the master SNP table
SUMMARY_MASTER_SNP_TABLE="out/summary_colocal/summary_master_SNP_table_grouped_by_piQTL.csv"
if [ ! -e ${SUMMARY_MASTER_SNP_TABLE} ]; then
    echo "Add colocaled pQTL, eQTL, and piQTL info into the master SNP table grouped by piQTL..."
    # Link to Jakobson et al pQTL data
    COLOCAL_PIQTL_EQTL_TABLE="out/piQTLs_vs_eQTLs/piQTLs_vs_eQTLs_colocal.csv"
    COLOCAL_PIQTL_PQTL_TABLE="out/piQTL_vs_pQTL/piQTLs_vs_pQTLs_colocal.csv"
    python src/summarize_colocaled_QTL_info_to_master_SNP_table_grouped_by_piQTL.py \
        --master_table ${MASTER_TABLE} \
        --piqtl_eqtl ${COLOCAL_PIQTL_EQTL_TABLE} \
        --piqtl_pqtl ${COLOCAL_PIQTL_PQTL_TABLE} \
        --output ${SUMMARY_MASTER_SNP_TABLE}
    echo "Done!"
fi


# Summarize colocaled piQTLs and eQTLs grouped by pQTLs into the master SNP table
SUMMARY_MASTER_SNP_TABLE_PQTL="out/summary_colocal/summary_master_SNP_table_grouped_by_pQTL.csv"
if [ ! -e ${SUMMARY_MASTER_SNP_TABLE_PQTL} ]; then
    echo "Add colocaled pQTL, eQTL, and piQTL info into the master SNP table grouped by pQTL..."
    # Link to Jakobson et al pQTL data
    COLOCAL_PIQTL_PQTL_TABLE="out/piQTL_vs_pQTL/piQTLs_vs_pQTLs_colocal.csv"
    COLOCAL_PQTL_EQTL_TABLE="out/pQTL_vs_eQTL/pQTLs_vs_eQTLs_colocal.csv"
    python src/summarize_colocaled_QTL_info_to_master_SNP_table_grouped_by_pQTL.py \
        --master_table ${MASTER_TABLE} \
        --piqtl_pqtl ${COLOCAL_PIQTL_PQTL_TABLE} \
        --pqtl_eqtl ${COLOCAL_PQTL_EQTL_TABLE} \
        --output ${SUMMARY_MASTER_SNP_TABLE_PQTL}
    echo "Done!"
fi


# Summarize colocaled piQTLs and pQTLs grouped by eQTLs into the master SNP table
SUMMARY_MASTER_SNP_TABLE_EQTL="out/summary_colocal/summary_master_SNP_table_grouped_by_eQTL.csv"
if [ ! -e ${SUMMARY_MASTER_SNP_TABLE_EQTL} ]; then
    echo "Add colocaled pQTL, eQTL, and piQTL info into the master SNP table grouped by eQTL..."
    # Link to Jakobson et al pQTL data
    COLOCAL_PIQTL_EQTL_TABLE="out/piQTLs_vs_eQTLs/piQTLs_vs_eQTLs_colocal.csv"
    COLOCAL_PQTL_EQTL_TABLE="out/pQTL_vs_eQTL/pQTLs_vs_eQTLs_colocal.csv"
    python src/summarize_colocaled_QTL_info_to_master_SNP_table_grouped_by_eQTL.py \
        --master_table ${MASTER_TABLE} \
        --piqtl_eqtl ${COLOCAL_PIQTL_EQTL_TABLE} \
        --pqtl_eqtl ${COLOCAL_PQTL_EQTL_TABLE} \
        --output ${SUMMARY_MASTER_SNP_TABLE_EQTL}
    echo "Done!"
fi
