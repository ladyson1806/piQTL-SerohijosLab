#!/bin/sh

cd $(dirname ${0})

mkdir -p out/formatted_tables

# Input files
TEYSSONNIERE_PQTL_INPUT="data/Teyssonniere_pQTLs.csv"
TEYSSONNIERE_EQTL_INPUT="data/Teyssonniere_eQTLs.csv"
PPI_GENE_LIST="data/PPI_gene_list.csv"



# Format pQTLs from Teyssonniere et al
TEYSSONNIERE_PQTL_FORMATTED="out/formatted_tables/Teyssonniere_pQTLs_formatted.csv"
if [ ! -e ${TEYSSONNIERE_PQTL_FORMATTED} ]; then
    echo "Formatting pQTLs from Teyssonniere et al..."
    python src/format_Teyssonniere_QTL_table.py \
        --input ${TEYSSONNIERE_PQTL_INPUT} \
        --output ${TEYSSONNIERE_PQTL_FORMATTED}
    echo "Done!"
fi

# Format eQTLs from Teyssonniere et al
TEYSSONNIERE_EQTL_FORMATTED="out/formatted_tables/Teyssonniere_eQTLs_formatted.csv"
if [ ! -e ${TEYSSONNIERE_EQTL_FORMATTED} ]; then
    echo "Formatting eQTLs from Teyssonniere et al..."
    python src/format_Teyssonniere_QTL_table.py \
        --input ${TEYSSONNIERE_EQTL_INPUT} \
        --output ${TEYSSONNIERE_EQTL_FORMATTED}
    echo "Done!"
fi


# Cerate a box plot of the cis and trans Teyssonniere pQTLs
TEYSSONNIERE_PQTL_BOXPLOT="out/figures/boxplot_Teyssonniere_pQTLs.png"
if [ ! -e ${TEYSSONNIERE_PQTL_BOXPLOT} ]; then
    echo "Create a box plot of the pQTLs from Teyssonniere et al..."
    python src/plot_boxplot.py \
        --input ${TEYSSONNIERE_PQTL_FORMATTED} \
        --output ${TEYSSONNIERE_PQTL_BOXPLOT} \
        --title "Teyssonniere pQTLs from natural isolates" \
        --y_label "abs(effect)" \
        --group_by "cis_trans" \
        --y_value "abs_effect"
    echo "Done!"
fi

## t-Test
TEYSSONNIERE_TTEST="out/t_test/Teyssonniere_pQTLs_t-test.txt"
if [ ! -e ${TEYSSONNIERE_TTEST} ]; then
    echo "Perform t-test on pQTLs (Teyssonniere et al)..."
    python src/perform_ttest.py \
        --input ${TEYSSONNIERE_PQTL_FORMATTED} \
        --output ${TEYSSONNIERE_TTEST} \
        --group_by "cis_trans" \
        --y_value "abs_effect" \
        --title "t-Test for the cis and trans Teyssonniere pQTLs"
    echo "Done!"
fi


# Create a box plot of the cis and trans Teyssonniere eQTLs
TEYSSONNIERE_EQTL_BOXPLOT="out/figures/boxplot_Teyssonniere_eQTLs.png"
if [ ! -e ${TEYSSONNIERE_EQTL_BOXPLOT} ]; then
    echo "Create a box plot of the eQTLs from Teyssonniere et al..."
    python src/plot_boxplot.py \
        --input ${TEYSSONNIERE_EQTL_FORMATTED} \
        --output ${TEYSSONNIERE_EQTL_BOXPLOT} \
        --title "Teyssonniere eQTLs from natural isolates" \
        --y_label "abs(effect)" \
        --group_by "cis_trans" \
        --y_value "abs_effect"
    echo "Done!"
fi

## t-Test
TEYSSONNIERE_EQTL_TTEST="out/t_test/Teyssonniere_eQTLs_t-test.txt"
if [ ! -e ${TEYSSONNIERE_EQTL_TTEST} ]; then
    echo "Perform t-test on eQTLs (Teyssonniere et al)..."
    python src/perform_ttest.py \
        --input ${TEYSSONNIERE_EQTL_FORMATTED} \
        --output ${TEYSSONNIERE_EQTL_TTEST} \
        --group_by "cis_trans" \
        --y_value "abs_effect" \
        --title "t-Test for the cis and trans Teyssonniere eQTLs"
    echo "Done!"
fi


##################################
# Analysis for our piQTL 44 genes
##################################

# pQTLs from Teyssonniere et al for our 44 piQTL genes
TEYSSONNIERE_PQTL_FORMATTED_PIQTL_GENES="out/formatted_tables/Teyssonniere_pQTLs_piQTL_genes.csv"
if [ ! -e ${TEYSSONNIERE_PQTL_FORMATTED_PIQTL_GENES} ]; then
    echo "Extract pQTLs for our 44 piQTL genes from Teyssonniere et al..."
    python src/format_Teyssonniere_QTL_table.py \
        --input ${TEYSSONNIERE_PQTL_INPUT} \
        --target_genes ${PPI_GENE_LIST} \
        --output ${TEYSSONNIERE_PQTL_FORMATTED_PIQTL_GENES}
    echo "Done!"
fi

## Cerate a box plot of the cis and trans pQTLs for our 44 piQTL genes
TEYSSONNIERE_PQTL_BOXPLOT_PIQTL_GENES="out/figures/boxplot_Teyssonniere_pQTLs_only_piQTL_genes.png"
if [ ! -e ${TEYSSONNIERE_PQTL_BOXPLOT_PIQTL_GENES} ]; then
    echo "Create a box plot of the pQTLs from Teyssonniere et al only for our 44 piQTL genes..."
    python src/plot_boxplot.py \
        --input ${TEYSSONNIERE_PQTL_FORMATTED_PIQTL_GENES} \
        --output ${TEYSSONNIERE_PQTL_BOXPLOT_PIQTL_GENES} \
        --title "pQTLs" \
        --y_label "abs(effect)" \
        --group_by "cis_trans" \
        --y_value "abs_effect"
    echo "Done!"
fi

## t-Test
TEYSSONNIERE_PQTL_TTEST_PIQTL_GENES="out/t_test/Teyssonniere_pQTLs_only_piQTL_genes_t-test.txt"
if [ ! -e ${TEYSSONNIERE_PQTL_TTEST_PIQTL_GENES} ]; then
    echo "Perform t-test on pQTLs (Teyssonniere et al)..."
    python src/perform_ttest.py \
        --input ${TEYSSONNIERE_PQTL_FORMATTED_PIQTL_GENES} \
        --output ${TEYSSONNIERE_PQTL_TTEST_PIQTL_GENES} \
        --group_by "cis_trans" \
        --y_value "abs_effect" \
        --title "t-Test for the cis and trans Teyssonniere pQTLs (only piQTL 44 genes)"
    echo "Done!"
fi


# eQTLs from Teyssonniere et al for our 44 piQTL genes
TEYSSONNIERE_EQTL_FORMATTED_PIQTL_GENES="out/formatted_tables/Teyssonniere_eQTLs_piQTL_genes.csv"
if [ ! -e ${TEYSSONNIERE_EQTL_FORMATTED_PIQTL_GENES} ]; then
    echo "Extract eQTLs for our 44 piQTL genes from Teyssonniere et al..."
    python src/format_Teyssonniere_QTL_table.py \
        --input ${TEYSSONNIERE_EQTL_INPUT} \
        --target_genes ${PPI_GENE_LIST} \
        --output ${TEYSSONNIERE_EQTL_FORMATTED_PIQTL_GENES}
    echo "Done!"
fi

## Cerate a box plot of the cis and trans eQTLs for our 44 piQTL genes
TEYSSONNIERE_EQTL_BOXPLOT_PIQTL_GENES="out/figures/boxplot_Teyssonniere_eQTLs_only_piQTL_genes.png"
if [ ! -e ${TEYSSONNIERE_EQTL_BOXPLOT_PIQTL_GENES} ]; then
    echo "Create a box plot of the eQTLs from Teyssonniere et al only for our 44 piQTL genes..."
    python src/plot_boxplot.py \
        --input ${TEYSSONNIERE_EQTL_FORMATTED_PIQTL_GENES} \
        --output ${TEYSSONNIERE_EQTL_BOXPLOT_PIQTL_GENES} \
        --title "eQTLs" \
        --y_label "abs(effect)" \
        --group_by "cis_trans" \
        --y_value "abs_effect"
    echo "Done!"
fi

## t-Test
TEYSSONNIERE_EQTL_TTEST_PIQTL_GENES="out/t_test/Teyssonniere_eQTLs_only_piQTL_genes_t-test.txt"
if [ ! -e ${TEYSSONNIERE_EQTL_TTEST_PIQTL_GENES} ]; then
    echo "Perform t-test on eQTLs (Teyssonniere et al)..."
    python src/perform_ttest.py \
        --input ${TEYSSONNIERE_EQTL_FORMATTED_PIQTL_GENES} \
        --output ${TEYSSONNIERE_EQTL_TTEST_PIQTL_GENES} \
        --group_by "cis_trans" \
        --y_value "abs_effect" \
        --title "t-Test for the cis and trans Teyssonniere eQTLs (only piQTL 44 genes)"
    echo "Done!"
fi
