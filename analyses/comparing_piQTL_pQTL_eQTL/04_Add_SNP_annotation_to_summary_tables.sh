#!/bin/sh

cd $(dirname ${0})

mkdir -p out/summary_annotated

# Link to SNP annotation table
SNP_ANNOTATION="data/SNP_annotation_table.csv"
if [ ! -e ${SNP_ANNOTATION} ]; then
    ln -s ../../../references/snp_metadata/snps_annotations_genome-version-3-64-1_with_gene_annotations_plus_LD_block_definition.txt ${SNP_ANNOTATION}
fi


# Merge SNP annotation into summary tables (Exact matched table)
SUMMARY_TABLE_EXACT="out/summary_annotated/summary_table_exact_matched_with_annotation.csv"
if [ ! -e ${SUMMARY_TABLE_EXACT} ]; then
    echo "Add SNP annotation into summary table of exact matched SNPs..."
    SUMMARY_TABLE_INPUT="out/summary/summary_master_SNP_tabl.csv"
    python src/add_SNP_annotation_to_summary_table.py \
        --summary_table ${SUMMARY_TABLE_INPUT} \
        --snp_annotation ${SNP_ANNOTATION} \
        --type exact \
        --output ${SUMMARY_TABLE_EXACT}
    echo "Done!"
fi


# Merge SNP annotation into summary tables (co-local matched table grouped by piQTL)
SUMMARY_TABLE_COLOCAL_PIQTL="out/summary_annotated/summary_table_colocaled_grouped_by_piQTL_with_annotation.csv"
if [ ! -e ${SUMMARY_TABLE_COLOCAL_PIQTL} ]; then
    echo "Add SNP annotation into summary table of colocaled SNPs grouped by piQTL..."
    SUMMARY_TABLE_INPUT="out/summary_colocal/summary_master_SNP_table_grouped_by_piQTL_full.csv"
    python src/add_SNP_annotation_to_summary_table.py \
        --summary_table ${SUMMARY_TABLE_INPUT} \
        --snp_annotation ${SNP_ANNOTATION} \
        --type colocal_groupby_piQTL \
        --output ${SUMMARY_TABLE_COLOCAL_PIQTL}
    echo "Done!"
fi


# Merge SNP annotation into summary tables (co-local matched table grouped by pQTL)
SUMMARY_TABLE_COLOCAL_PQTL="out/summary_annotated/summary_table_colocaled_grouped_by_pQTL_with_annotation.csv"
if [ ! -e ${SUMMARY_TABLE_COLOCAL_PQTL} ]; then
    echo "Add SNP annotation into summary table of colocaled SNPs grouped by pQTL..."
    SUMMARY_TABLE_INPUT="out/summary_colocal/summary_master_SNP_table_grouped_by_pQTL_full.csv"
    python src/add_SNP_annotation_to_summary_table.py \
        --summary_table ${SUMMARY_TABLE_INPUT} \
        --snp_annotation ${SNP_ANNOTATION} \
        --type colocal_groupby_pQTL \
        --output ${SUMMARY_TABLE_COLOCAL_PQTL}
    echo "Done!"
fi


# Merge SNP annotation into summary tables (co-local matched table grouped by eQTL)
SUMMARY_TABLE_COLOCAL_EQTL="out/summary_annotated/summary_table_colocaled_grouped_by_eQTL_with_annotation.csv"
if [ ! -e ${SUMMARY_TABLE_COLOCAL_EQTL} ]; then
    echo "Add SNP annotation into summary table of colocaled SNPs grouped by eQTL..."
    SUMMARY_TABLE_INPUT="out/summary_colocal/summary_master_SNP_table_grouped_by_eQTL_full.csv"
    python src/add_SNP_annotation_to_summary_table.py \
        --summary_table ${SUMMARY_TABLE_INPUT} \
        --snp_annotation ${SNP_ANNOTATION} \
        --type colocal_groupby_eQTL \
        --output ${SUMMARY_TABLE_COLOCAL_EQTL}
    echo "Done!"
fi
