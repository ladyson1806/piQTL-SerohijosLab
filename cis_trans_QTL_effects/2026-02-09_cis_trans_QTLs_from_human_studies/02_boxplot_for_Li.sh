#!/bin/sh

cd $(dirname ${0})

mkdir -p out/formatted_tables

# Input files
LI_PQTL_INPUT="data/Li_pQTLs.csv"


# Format pQTLs from Li et al
LI_PQTL_FORMATTED="out/formatted_tables/Li_pQTLs_formatted.csv"
if [ ! -e ${LI_PQTL_FORMATTED} ]; then
    echo "Formatting pQTLs from Li et al..."
    python src/format_Li_QTL_table.py \
        --input ${LI_PQTL_INPUT} \
        --output ${LI_PQTL_FORMATTED}
    echo "Done!"
fi

# Cerate a box plot of the cis and trans Li pQTLs
LI_PQTL_BOXPLOT="out/figures/boxplot_Li_pQTLs.png"
if [ ! -e ${LI_PQTL_BOXPLOT} ]; then
    echo "Create a box plot of the pQTLs from Li et al..."
    python src/plot_boxplot.py \
        --input ${LI_PQTL_FORMATTED} \
        --output ${LI_PQTL_BOXPLOT} \
        --title "Li pQTLs from natural isolates" \
        --y_label "abs(effect)" \
        --group_by "cis_trans" \
        --y_value "abs_effect"
    echo "Done!"
fi

## t-Test
LI_TTEST="out/t_test/Li_pQTLs_t-test.txt"
if [ ! -e ${LI_TTEST} ]; then
    echo "Perform t-test on pQTLs (Li et al)..."
    python src/perform_ttest.py \
        --input ${LI_PQTL_FORMATTED} \
        --output ${LI_TTEST} \
        --group_by "cis_trans" \
        --y_value "abs_effect" \
        --title "t-Test for the cis and trans Li pQTLs"
    echo "Done!"
fi
