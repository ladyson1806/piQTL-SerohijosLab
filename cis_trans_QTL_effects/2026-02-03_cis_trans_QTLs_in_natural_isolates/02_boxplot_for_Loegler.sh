#!/bin/sh

cd $(dirname ${0})

mkdir -p out/formatted_tables

# Input files
LOEGLER_QTL_INPUT="data/Loegler_QTLs.csv" # Loegler et al make a combined table of eQTLs and pQTLs, and growth QTLs
PPI_GENE_LIST="data/PPI_gene_list.csv"



# Format pQTLs from Loegler et al
LOEGLER_PQTL_FORMATTED="out/formatted_tables/Loegler_pQTLs_formatted.csv"
if [ ! -e ${LOEGLER_PQTL_FORMATTED} ]; then
    echo "Formatting pQTLs from Loegler et al..."
    python src/format_Loegler.py \
        --input ${LOEGLER_QTL_INPUT} \
        --pheno_type "Proteomics" \
        --output ${LOEGLER_PQTL_FORMATTED}
    echo "Done!"
fi


# Format eQTLs from Loegler et al
LOEGLER_EQTL_FORMATTED="out/formatted_tables/Loegler_eQTLs_formatted.csv"
if [ ! -e ${LOEGLER_EQTL_FORMATTED} ]; then
    echo "Formatting eQTLs from Loegler et al..."
    python src/format_Loegler.py \
        --input ${LOEGLER_QTL_INPUT} \
        --pheno_type "Transcriptomics" \
        --output ${LOEGLER_EQTL_FORMATTED}
    echo "Done!"
fi


# Cerate a box plot of the cis and trans Loegler pQTLs
LOEGLER_PQTL_BOXPLOT="out/figures/boxplot_Loegler_pQTLs.png"
if [ ! -e ${LOEGLER_PQTL_BOXPLOT} ]; then
    echo "Create a box plot of the pQTLs from Loegler et al..."
    python src/plot_boxplot.py \
        --input ${LOEGLER_PQTL_FORMATTED} \
        --output ${LOEGLER_PQTL_BOXPLOT} \
        --title "Loegler pQTLs from natural isolates" \
        --y_label "abs(effect)" \
        --group_by "cis_trans" \
        --y_value "abs_effect"
    echo "Done!"
fi

## t-Test
LOEGLER_TTEST="out/t_test/Loegler_pQTLs_t-test.txt"
if [ ! -e ${LOEGLER_TTEST} ]; then
    echo "Perform t-test on pQTLs (Loegler et al)..."
    python src/perform_ttest.py \
        --input ${LOEGLER_PQTL_FORMATTED} \
        --output ${LOEGLER_TTEST} \
        --group_by "cis_trans" \
        --y_value "abs_effect" \
        --title "t-Test for the cis and trans Loegler pQTLs"
    echo "Done!"
fi


# Create a box plot of the cis and trans Loegler eQTLs
LOEGLER_EQTL_BOXPLOT="out/figures/boxplot_Loegler_eQTLs.png"
if [ ! -e ${LOEGLER_EQTL_BOXPLOT} ]; then
    echo "Create a box plot of the eQTLs from Loegler et al..."
    python src/plot_boxplot.py \
        --input ${LOEGLER_EQTL_FORMATTED} \
        --output ${LOEGLER_EQTL_BOXPLOT} \
        --title "Loegler eQTLs from natural isolates" \
        --y_label "abs(effect)" \
        --group_by "cis_trans" \
        --y_value "abs_effect"
    echo "Done!"
fi

## t-Test
LOEGLER_EQTL_TTEST="out/t_test/Loegler_eQTLs_t-test.txt"
if [ ! -e ${LOEGLER_EQTL_TTEST} ]; then
    echo "Perform t-test on eQTLs (Loegler et al)..."
    python src/perform_ttest.py \
        --input ${LOEGLER_EQTL_FORMATTED} \
        --output ${LOEGLER_EQTL_TTEST} \
        --group_by "cis_trans" \
        --y_value "abs_effect" \
        --title "t-Test for the cis and trans Loegler eQTLs"
    echo "Done!"
fi


##################################
# Analysis for our piQTL 44 genes
##################################

# pQTLs from Loegler et al for our 44 piQTL genes
LOEGLER_PQTL_FORMATTED_PIQTL_GENES="out/formatted_tables/Loegler_pQTLs_piQTL_genes.csv"
if [ ! -e ${LOEGLER_PQTL_FORMATTED_PIQTL_GENES} ]; then
    echo "Extract pQTLs for our 44 piQTL genes from Loegler et al..."
    python src/format_Loegler.py \
        --input ${LOEGLER_QTL_INPUT} \
        --pheno_type "Proteomics" \
        --target_genes ${PPI_GENE_LIST} \
        --output ${LOEGLER_PQTL_FORMATTED_PIQTL_GENES}
    echo "Done!"
fi

## Cerate a box plot of the cis and trans pQTLs for our 44 piQTL genes
LOEGLER_PQTL_BOXPLOT_PIQTL_GENES="out/figures/boxplot_Loegler_pQTLs_only_piQTL_genes.png"
if [ ! -e ${LOEGLER_PQTL_BOXPLOT_PIQTL_GENES} ]; then
    echo "Create a box plot of the pQTLs from Loegler et al only for our 44 piQTL genes..."
    python src/plot_boxplot.py \
        --input ${LOEGLER_PQTL_FORMATTED_PIQTL_GENES} \
        --output ${LOEGLER_PQTL_BOXPLOT_PIQTL_GENES} \
        --title "pQTLs" \
        --y_label "abs(effect)" \
        --group_by "cis_trans" \
        --y_value "abs_effect"
    echo "Done!"
fi

## t-Test
LOEGLER_PQTL_TTEST_PIQTL_GENES="out/t_test/Loegler_pQTLs_only_piQTL_genes_t-test.txt"
if [ ! -e ${LOEGLER_PQTL_TTEST_PIQTL_GENES} ]; then
    echo "Perform t-test on pQTLs (Loegler et al)..."
    python src/perform_ttest.py \
        --input ${LOEGLER_PQTL_FORMATTED_PIQTL_GENES} \
        --output ${LOEGLER_PQTL_TTEST_PIQTL_GENES} \
        --group_by "cis_trans" \
        --y_value "abs_effect" \
        --title "t-Test for the cis and trans Loegler pQTLs (only piQTL 44 genes)"
    echo "Done!"
fi


# eQTLs from Loegler et al for our 44 piQTL genes
LOEGLER_EQTL_FORMATTED_PIQTL_GENES="out/formatted_tables/Loegler_eQTLs_piQTL_genes.csv"
if [ ! -e ${LOEGLER_EQTL_FORMATTED_PIQTL_GENES} ]; then
    echo "Extract eQTLs for our 44 piQTL genes from Loegler et al..."
    python src/format_Loegler.py \
        --input ${LOEGLER_QTL_INPUT} \
        --pheno_type "Transcriptomics" \
        --target_genes ${PPI_GENE_LIST} \
        --output ${LOEGLER_EQTL_FORMATTED_PIQTL_GENES}
    echo "Done!"
fi

## Cerate a box plot of the cis and trans eQTLs for our 44 piQTL genes
LOEGLER_EQTL_BOXPLOT_PIQTL_GENES="out/figures/boxplot_Loegler_eQTLs_only_piQTL_genes.png"
if [ ! -e ${LOEGLER_EQTL_BOXPLOT_PIQTL_GENES} ]; then
    echo "Create a box plot of the eQTLs from Loegler et al only for our 44 piQTL genes..."
    python src/plot_boxplot.py \
        --input ${LOEGLER_EQTL_FORMATTED_PIQTL_GENES} \
        --output ${LOEGLER_EQTL_BOXPLOT_PIQTL_GENES} \
        --title "eQTLs" \
        --y_label "abs(effect)" \
        --group_by "cis_trans" \
        --y_value "abs_effect"
    echo "Done!"
fi

## t-Test
LOEGLER_EQTL_TTEST_PIQTL_GENES="out/t_test/Loegler_eQTLs_only_piQTL_genes_t-test.txt"
if [ ! -e ${LOEGLER_EQTL_TTEST_PIQTL_GENES} ]; then
    echo "Perform t-test on eQTLs (Loegler et al)..."
    python src/perform_ttest.py \
        --input ${LOEGLER_EQTL_FORMATTED_PIQTL_GENES} \
        --output ${LOEGLER_EQTL_TTEST_PIQTL_GENES} \
        --group_by "cis_trans" \
        --y_value "abs_effect" \
        --title "t-Test for the cis and trans Loegler eQTLs (only piQTL 44 genes)"
    echo "Done!"
fi
