#!/bin/sh

cd $(dirname ${0})

mkdir -p out/formatted_tables

# Input files
VOSA_CIS_EQTL_INPUT="data/Vosa_cis_eQTLs.tsv"
VOSA_TRANS_EQTL_INPUT="data/Vosa_trans_eQTLs.tsv"


# Format pQTLs from Li et al
VOSA_EQTL_FORMATTED="out/formatted_tables/Vosa_eQTLs_formatted.csv"
if [ ! -e ${VOSA_EQTL_FORMATTED} ]; then
    echo "Formatting eQTLs from Vosa et al..."
    python src/format_Vosa_QTL_table.py \
        --input_cis ${VOSA_CIS_EQTL_INPUT} \
        --input_trans ${VOSA_TRANS_EQTL_INPUT} \
        --output ${VOSA_EQTL_FORMATTED}
    echo "Done!"
fi


# Cerate a box plot of the cis and trans Vosa eQTLs
VOSA_EQTL_BOXPLOT="out/figures/boxplot_Vosa_eQTLs.png"
if [ ! -e ${VOSA_EQTL_BOXPLOT} ]; then
    echo "Create a box plot of the eQTLs from Vosa et al..."
    python src/plot_boxplot.py \
        --input ${VOSA_EQTL_FORMATTED} \
        --output ${VOSA_EQTL_BOXPLOT} \
        --title "Vosa eQTLs from human blood samples" \
        --y_label "|Z|" \
        --group_by "cis_trans" \
        --y_value "abs_z"
    echo "Done!"
fi

# Make a zoomed-in version of the boxplot
VOSA_EQTL_BOXPLOT_ZOOMED="out/figures/boxplot_Vosa_eQTLs_zoomed.png"
if [ ! -e ${VOSA_EQTL_BOXPLOT_ZOOMED} ]; then
    echo "Create a zoomed-in box plot of the eQTLs from Vosa et al..."
    python src/plot_boxplot.py \
        --input ${VOSA_EQTL_FORMATTED} \
        --output ${VOSA_EQTL_BOXPLOT_ZOOMED} \
        --title "Vosa eQTLs from human blood samples (zoomed-in)" \
        --y_label "|Z|" \
        --group_by "cis_trans" \
        --y_value "abs_z" \
        --y_limit 80
    echo "Done!"
fi


## t-Test
VOSA_TTEST="out/t_test/Vosa_eQTLs_t-test.txt"
if [ ! -e ${VOSA_TTEST} ]; then
    echo "Perform t-test on eQTLs (Vosa et al)..."
    python src/perform_ttest.py \
        --input ${VOSA_EQTL_FORMATTED} \
        --output ${VOSA_TTEST} \
        --group_by "cis_trans" \
        --y_value "abs_z" \
        --title "t-Test for the cis and trans Vosa eQTLs"
    echo "Done!"
fi
