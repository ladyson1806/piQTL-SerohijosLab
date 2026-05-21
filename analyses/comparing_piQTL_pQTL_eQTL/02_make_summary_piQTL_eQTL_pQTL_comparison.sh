#!/bin/sh

cd $(dirname ${0})

mkdir -p out/summary

# Link to Albert's genotype table
ALBERT_GENOTYPE="data/Albert_genotype_table.csv"
if [ ! -e ${ALBERT_GENOTYPE} ]; then
    ln -s ../../../references/Albert_eLife_2018/SI_Data_03_genotypes_header_tab_to_enter.txt ${ALBERT_GENOTYPE}
fi

# Prepare template SNP table for storing pi/e/pQTLs.
COMMON_SNPS_PIQTL_ALBERT="out/summary/master_SNP_table.csv"
if [ ! -e ${COMMON_SNPS_PIQTL_ALBERT} ]; then
    echo "Construct master SNP table summarizing SNPs used in both of Albert and Besse..."
    python src/construct_combined_SNP_table.py \
        --piqtl data/Jakobson_ld_info.csv \
        --eqtl data/Albert_genotype_table.csv \
        --output ${COMMON_SNPS_PIQTL_ALBERT}
    echo "Done!"
fi


# Add pQTL, eQTL, and piQTL info into the master SNP table
SUMMARY_MASTER_SNP_TABLE="out/summary/summary_master_SNP_tabl.csv"
if [ ! -e ${SUMMARY_MASTER_SNP_TABLE} ]; then
    echo "Add pQTL, eQTL, and piQTL info into the master SNP table..."
    # Link to Jakobson et al pQTL data
    PIQTL_TABLE="data/piQTLs_formatted_lead.csv"
    EQTL_TABLE="data/eQTLs_formatted.csv"
    PQTL_TABLE="out/formatted_tables/Jakobson_pQTLs_with_LD_blocks.csv"
    python src/add_QTL_info_to_master_SNP_table.py \
        --master_table ${COMMON_SNPS_PIQTL_ALBERT} \
        --piqtl ${PIQTL_TABLE} \
        --eqtl ${EQTL_TABLE} \
        --pqtl ${PQTL_TABLE} \
        --output ${SUMMARY_MASTER_SNP_TABLE}
    echo "Done!"
fi

exit 0

SUMMARY_UNIQUE_SNPS="out/summary/unique_snps_piQTL_pQTL_eQTL.txt"
if [ ! -e ${SUMMARY_UNIQUE_SNPS} ]; then
    echo "Extract unique SNPs that are piQTL, pQTL, and eQTL..."
    python src/extract_unique_snps_piQTL_pQTL_eQTL.py \
        --piqtl data/piQTLs_formatted_lead.csv \
        --pqtl ${JAKOBSON_PQTL_WITH_LD} \
        --eqtl data/eQTLs_formatted.csv \
        --output ${SUMMARY_UNIQUE_SNPS}

    echo "Done!"
fi

## Total SNPs for piQTLs and pQTLs (Jakobson et al) are 12054 SNPs

# Extract unique SNPs ID from piQTLs


## From there,
