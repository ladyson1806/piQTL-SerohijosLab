#!/bin/sh

cd $(dirname ${0})

mkdir -p data
mkdir -p out/formatted_tables
mkdir -p out/figures
mkdir -p out/t_test

# Preprocessing: Reformat Jakobson's data for boxplot and t-test
JAKOBSON_INPUT="data/Jakobson_pQTLs.csv"

# Identify cis/trans-pQTLs from Jakobson et al
JAKOBSON_CIS_TRANS="out/formatted_tables/Jakobson_pQTLs.csv"
if [ ! -e ${JAKOBSON_CIS_TRANS} ]; then
    echo "Identify cis- and trans-pQTLs from Jakobson et al..."
    python src/identify_cis_pQTLs_from_Jakobson.py \
        --input ${JAKOBSON_INPUT} \
        --output ${JAKOBSON_CIS_TRANS}
    echo "Done!"
fi

# Cerate a box plot of the cis and trans pQTLs
JAKOBSON_BOXPLOT="out/figures/boxplot_Jakobson_pQTLs.png"
if [ ! -e ${JAKOBSON_BOXPLOT} ]; then
    echo "Create a box plot of the pQTLs from Jakobson et al..."
    python src/plot_boxplot.py \
        --input ${JAKOBSON_CIS_TRANS} \
        --output ${JAKOBSON_BOXPLOT} \
        --title "pQTLs" \
        --y_label "abs(beta)" \
        --group_by "cis_trans" \
        --y_value "abs_beta"
    echo "Done!"
fi

# t-Test
JAKOBSON_TTEST="out/t_test/Jakobson_pQTLs_t-test.txt"
if [ ! -e ${JAKOBSON_TTEST} ]; then
    echo "Perform t-test on pQTLs (Jakobson et al)..."
    python src/perform_ttest.py \
        --input ${JAKOBSON_CIS_TRANS} \
        --output ${JAKOBSON_TTEST} \
        --group_by "cis_trans" \
        --y_value "abs_beta" \
        --title "t-Test for the cis and trans Jakobson pQTLs"
    echo "Done!"
fi


##################################
# Analysis for our piQTL 44 genes
##################################
JAKOBSON_CIS_TRANS_PIQTL_GENES="out/formatted_tables/Jakobson_pQTLs_piQTL_genes.csv"
if [ ! -e ${JAKOBSON_CIS_TRANS_PIQTL_GENES} ]; then
    echo "Extract pQTLs for our 44 piQTL genes from Jakobson et al..."
    TARGET_GENES="data/2025-04-17_PPI_gene_list.csv"
    python src/identify_cis_pQTLs_from_Jakobson.py \
        --input ${JAKOBSON_INPUT} \
        --target_genes ${TARGET_GENES} \
        --output ${JAKOBSON_CIS_TRANS_PIQTL_GENES}
    echo "Done!"
fi

# Cerate a box plot of the cis and trans pQTLs
JAKOBSON_BOXPLOT_PIQTL_GENE="out/figures/boxplot_Jakobson_pQTLs_only_piQTL_genes.png"
if [ ! -e ${JAKOBSON_BOXPLOT_PIQTL_GENE} ]; then
    echo "Create a box plot of the pQTLs from Jakobson et al only for our 44 piQTL genes..."
    python src/plot_boxplot.py \
        --input ${JAKOBSON_CIS_TRANS_PIQTL_GENES} \
        --output ${JAKOBSON_BOXPLOT_PIQTL_GENE} \
        --title "pQTLs" \
        --y_label "abs(beta)" \
        --group_by "cis_trans" \
        --y_value "abs_beta"
    echo "Done!"
fi

# t-Test
JAKOBSON_TTEST_PIQTL_GENE="out/t_test/Jakobson_pQTLs_only_piQTL_genes_t-test.txt"
if [ ! -e ${JAKOBSON_TTEST_PIQTL_GENE} ]; then
    echo "Perform t-test on pQTLs (Jakobson et al)..."
    python src/perform_ttest.py \
        --input ${JAKOBSON_CIS_TRANS_PIQTL_GENES} \
        --output ${JAKOBSON_TTEST_PIQTL_GENE} \
        --group_by "cis_trans" \
        --y_value "abs_beta" \
        --title "t-Test for the cis and trans Jakobson pQTLs (only piQTL 44 genes)"
    echo "Done!"
fi
