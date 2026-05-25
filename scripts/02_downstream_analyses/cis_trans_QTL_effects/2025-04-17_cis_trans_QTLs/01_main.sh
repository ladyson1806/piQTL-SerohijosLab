#!/bin/sh

cd $(dirname ${0})


# Identify cis/trans-piQTLs
OUTPUT="out/formatted_tables/piQTLs.csv"
if [ ! -e ${OUTPUT} ]; then
    echo "Identify cis- and trans-piQTLs..."
    INPUT="data/2025-04-17_piQTL_SNPs.csv"
    ANNOTATIONS="data/2025-04-17_snps_annotations_LD050_cisSNP-5kb.csv"
    python src/identify_cis_piQTLs.py \
        --input ${INPUT} \
        --annotations ${ANNOTATIONS} \
        --output ${OUTPUT}
    echo "Done!"
fi

# Identify cis/trans-eQTLs from Albert et al only focusing on our PPI genes
OUTPUT="out/formatted_tables/Albert_eQTLs.csv"
if [ ! -e ${OUTPUT} ]; then
    echo "Identify cis- and trans-eQTLs from Albert et al..."
    TARGET_GENES="data/2025-04-17_PPI_gene_list.csv"
    python src/identify_cis_eQTLs_from_Albert.py \
        --input ${INPUT} \
        --target_genes ${TARGET_GENES} \
        --output ${OUTPUT}
    echo "Done!"
fi


# Identify cis/trans pQTL from Chick et al
INPUT="data/2025-04-17_Chick_pQTL_Table_S5.csv"
OUTPUT="out/formatted_tables/Chick_pQTLs.csv"
if [ ! -e ${OUTPUT} ]; then
    echo "Identify cis- and trans-pQTLs from Chick et al..."
    python src/identify_cis_pQTLs_from_Chick.py \
        --input ${INPUT} \
        --output ${OUTPUT}
    echo "Done!"
fi


# Identify cis/trans eQTL from Chick et al
INPUT="data/2025-04-18_Chick_eQTL_Table_S4.csv"
OUTPUT="out/formatted_tables/Chick_eQTLs.csv"
if [ ! -e ${OUTPUT} ]; then
    echo "Identify cis- and trans-eQTLs from Chick et al..."
    python src/identify_cis_eQTLs_from_Chick.py \
        --input ${INPUT} \
        --output ${OUTPUT}
    echo "Done!"
fi


##################
# Box plot
##################

# Cerate a box plot of the cis and trans QTLs
## piQTLs
INPUT="out/formatted_tables/piQTLs.csv"
OUTPUT="out/figures/piQTLs_boxplot.png"
if [ ! -e ${OUTPUT} ]; then echo "Create a box plot of the piQTLs..."
    python src/plot_boxplot.py \
        --input ${INPUT} \
        --output ${OUTPUT} \
        --title "piQTLs" \
        --y_label "abs(beta)" \
        --group_by "cis_trans" \
        --y_value "abs_beta"
    echo "Done!"
fi


## eQTLs
INPUT="out/formatted_tables/Albert_eQTLs.csv"
OUTPUT="out/figures/Albert_eQTLs_boxplot.png"

if [ ! -e ${OUTPUT} ]; then echo "Create a box plot of the eQTLs..."
    python src/plot_boxplot.py \
        --input ${INPUT} \
        --output ${OUTPUT} \
        --title "eQTLs" \
        --y_label "abs(beta)" \
        --group_by "cis_trans" \
        --y_value "beta"
    echo "Done!"
fi


## eQTLs (all genes)
INPUT="out/formatted_tables/Albert_eQTLs_all.csv"
OUTPUT="out/figures/Albert_eQTLs_all_boxplot.png"
if [ ! -e ${OUTPUT} ]; then echo "Create a box plot of the eQTLs (all genes)..."
    python src/plot_boxplot.py \
        --input ${INPUT} \
        --output ${OUTPUT} \
        --title "eQTLs (all genes)" \
        --y_label "abs(beta)" \
        --group_by "cis_trans" \
        --y_value "beta"
    echo "Done!"
fi


## pQTLs from Chick et al
INPUT="out/formatted_tables/Chick_pQTLs.csv"
OUTPUT="out/figures/Chick_pQTLs_boxplot.png"
if [ ! -e ${OUTPUT} ]; then echo "Create a box plot of the pQTLs..."
    python src/plot_boxplot.py \
        --input ${INPUT} \
        --output ${OUTPUT} \
        --title "pQTLs (Chick et al; Inbred mouse pQTL)" \
        --y_label "abs(beta)" \
        --group_by "cis_trans" \
        --y_value "beta"
    echo "Done!"
fi


## eQTLs from Chick et al
INPUT="out/formatted_tables/Chick_eQTLs.csv"
OUTPUT="out/figures/Chick_eQTLs_boxplot.png"
if [ ! -e ${OUTPUT} ]; then echo "Create a box plot of the eQTLs (Chick et al)..."
    python src/plot_boxplot.py \
        --input ${INPUT} \
        --output ${OUTPUT} \
        --title "eQTLs (Chick et al; Inbred mouse eQTL)" \
        --y_label "abs(beta)" \
        --group_by "cis_trans" \
        --y_value "beta"
    echo "Done!"
fi


#####
# t-tests
####

# piQTLs
INPUT="out/formatted_tables/piQTLs.csv"
OUTPUT="out/t_test/piQTLs_t-test.txt"
if [ ! -e ${OUTPUT} ]; then echo "Perform t-test on piQTLs..."
    python src/perform_ttest.py \
        --input ${INPUT} \
        --output ${OUTPUT} \
        --group_by "cis_trans" \
        --y_value "abs_beta" \
        --title "t-Test for the cis and trans piQTLs"
    echo "Done!"
fi


# eQTLs
INPUT="out/formatted_tables/Albert_eQTLs.csv"
OUTPUT="out/t_test/Albert_eQTLs_t-test.txt"
if [ ! -e ${OUTPUT} ]; then echo "Perform t-test on eQTLs..."
    python src/perform_ttest.py \
        --input ${INPUT} \
        --output ${OUTPUT} \
        --group_by "cis_trans" \
        --y_value "beta" \
        --title "t-Test for the cis and trans Albert eQTLs (only PPI-geenes)"
    echo "Done!"
fi

# eQTLs (all genes)
INPUT="out/formatted_tables/Albert_eQTLs_all.csv"
OUTPUT="out/t_test/Albert_eQTLs_all_t-test.txt"
if [ ! -e ${OUTPUT} ]; then echo "Perform t-test on eQTLs (all genes)..."
    python src/perform_ttest.py \
        --input ${INPUT} \
        --output ${OUTPUT} \
        --group_by "cis_trans" \
        --y_value "beta" \
        --title "t-Test for the cis and trans Albert eQTLs (all genes)"
    echo "Done!"
fi


# pQTLs from Chick et al
INPUT="out/formatted_tables/Chick_pQTLs.csv"
OUTPUT="out/t_test/Chick_pQTLs_t-test.txt"
if [ ! -e ${OUTPUT} ]; then echo "Perform t-test on pQTLs..."
    python src/perform_ttest.py \
        --input ${INPUT} \
        --output ${OUTPUT} \
        --group_by "cis_trans" \
        --y_value "beta" \
        --title "t-Test for the cis and trans Chick pQTLs (Inbred mouse pQTL)"
    echo "Done!"
fi


# eQTLs from Chick et al
INPUT="out/formatted_tables/Chick_eQTLs.csv"
OUTPUT="out/t_test/Chick_eQTLs_t-test.txt"
if [ ! -e ${OUTPUT} ]; then echo "Perform t-test on eQTLs (Chick et al)..."
    python src/perform_ttest.py \
        --input ${INPUT} \
        --output ${OUTPUT} \
        --group_by "cis_trans" \
        --y_value "beta" \
        --title "t-Test for the cis and trans Chick eQTLs (Inbred mouse eQTL)"
    echo "Done!"
fi
