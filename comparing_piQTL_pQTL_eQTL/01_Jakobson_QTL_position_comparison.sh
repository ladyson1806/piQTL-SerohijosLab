#!/bin/sh

cd $(dirname ${0})

mkdir -p data
mkdir -p out/formatted_tables
mkdir -p out/figures
mkdir -p out/t_test

# Preprocessing: Parsing LB block range of Jakobson's data
JAKOBSON_LD_RANGE="out/formatted_tables/Jakobson_LD_block_ranges.csv"
if [ ! -e ${JAKOBSON_LD_RANGE} ]; then
    echo "Parse LB block ranges from Jakobson et al..."

    # Define paths to Jakobson et al data
    ## LD table
    JAKOBSON_LD_TABLE="data/Jakobson_ld_table.txt"
    ## LD info
    JAKOBSON_LD_INFO="data/Jakobson_ld_info.csv"

    # Parse LD block ranges
    python src/parse_Jakobson_LD_block_ranges.py \
        --ld_table ${JAKOBSON_LD_TABLE} \
        --ld_info ${JAKOBSON_LD_INFO} \
        --output ${JAKOBSON_LD_RANGE}
    echo "Done!"
fi


# Add LD block range info into Jakobson pQTL data
JAKOBSON_PQTL_WITH_LD="out/formatted_tables/Jakobson_pQTLs_with_LD_blocks.csv"
if [ ! -e ${JAKOBSON_PQTL_WITH_LD} ]; then
    echo "Add LD block range info into Jakobson pQTL data..."
    # Link to Jakobson et al pQTL data
    JAKOBSON_PQTL="out/formatted_tables/Jakobson_pQTLs.csv"

    python src/add_LD_block_info_to_Jakobson_pQTLs.py \
        --pQTL_input ${JAKOBSON_PQTL} \
        --ld_range_input ${JAKOBSON_LD_RANGE} \
        --output ${JAKOBSON_PQTL_WITH_LD}
    echo "Done!"
fi


# Comparing Jakobson pQTLs with our piQTLs
JAKOBSON_PQTL_PIQTL_OVERLAP="out/piQTL_vs_pQTL/piQTLs_vs_pQTLs_exact.csv"
if [ ! -e ${JAKOBSON_PQTL_PIQTL_OVERLAP} ]; then
    echo "Extract overlapped pQTLs (Jakobson et al.) with our piQTLs..."
    python src/extract_overlapped_pQTLs_with_piQTLs.py \
        --piqtl data/piQTLs_formatted_lead.csv \
        --pqtl ${JAKOBSON_PQTL_WITH_LD} \
        --output_dir out/piQTL_vs_pQTL
    echo "Done!"
fi


# Comparing Jakobson pQTLs with Albert eQTLs
JAKOBSON_PQTL_EQTL_OVERLAP="out/pQTL_vs_eQTL/pQTLs_vs_eQTLs_exact.csv"
if [ ! -e ${JAKOBSON_PQTL_EQTL_OVERLAP} ]; then
    echo "Extract overlapped pQTLs (Jakobson et al.) with Albert eQTLs..."
    python src/extract_overlapped_pQTLs_with_eQTLs.py \
        --eqtl data/eQTLs_formatted.csv \
        --pqtl ${JAKOBSON_PQTL_WITH_LD} \
        --output_dir out/pQTL_vs_eQTL
    echo "Done!"
fi


# Comparing Jakobson pQTLs with piQTL-eQTL already overlapped SNPs
# In other words, identify three-QTL overlapping SNPs
JAKOBSON_PQTL_PIQTL_EQTL_OVERLAP="out/piQTL_vs_eQTL_vs_pQTL/piQTLs_vs_eQTLs_vs_pQTLs_exact.csv"
if [ ! -e ${JAKOBSON_PQTL_PIQTL_EQTL_OVERLAP} ]; then
    echo "Extract overlapped pQTLs (Jakobson et al.) with piQTL-eQTL and eQTL overlapped SNPs..."
    python src/extract_overlapped_pQTLs_with_piQTLs_eQTLs.py \
        --piqtl_eqtl out/piQTLs_vs_eQTLs/piQTLs_vs_eQTLs_colocal.csv \
        --pqtl ${JAKOBSON_PQTL_WITH_LD} \
        --output_dir out/piQTL_vs_eQTL_vs_pQTL
    echo "Done!"
fi



# Extract pQTLs that directly affect piQTL genes
## Exact matched cases
JAKOBSON_PQTL_PIQTL_GENE_TARGETED="out/piQTL_vs_pQTL/piQTLs_directly_affected_exact.csv"
if [ ! -e ${JAKOBSON_PQTL_PIQTL_GENE_TARGETED} ]; then
    echo "Extract pQTLs (Jakobson et al.) that directly affect piQTL genes..."
    python src/extract_directly_affected_pQTLs.py \
        --data ${JAKOBSON_PQTL_PIQTL_OVERLAP} \
        --output ${JAKOBSON_PQTL_PIQTL_GENE_TARGETED}
    echo "Done!"
fi

## Co-localized cases
JAKOBSON_PQTL_PIQTL_GENE_TARGETED_COLOC="out/piQTL_vs_pQTL/piQTLs_directly_affected_colocal.csv"
if [ ! -e ${JAKOBSON_PQTL_PIQTL_GENE_TARGETED_COLOC} ]; then
    echo "Extract pQTLs (Jakobson et al.) that directly affect piQTL genes (colocalized)..."
    python src/extract_directly_affected_pQTLs.py \
        --data out/piQTL_vs_pQTL/piQTLs_vs_pQTLs_colocal.csv \
        --output ${JAKOBSON_PQTL_PIQTL_GENE_TARGETED_COLOC}
    echo "Done!"
fi


# Extract pQTLs that directly affect eQTL genes
## Exact matched cases
JAKOBSON_PQTL_EQTL_GENE_TARGETED="out/pQTL_vs_eQTL/pQTLs_directly_affected_exact.csv"
if [ ! -e ${JAKOBSON_PQTL_EQTL_GENE_TARGETED} ]; then
    echo "Extract pQTLs (Jakobson et al.) that directly affect eQTL genes..."
    python src/extract_directly_affected_eQTLs.py \
        --data ${JAKOBSON_PQTL_EQTL_OVERLAP} \
        --output ${JAKOBSON_PQTL_EQTL_GENE_TARGETED}
    echo "Done!"
fi

## Co-localized cases
JAKOBSON_PQTL_EQTL_GENE_TARGETED_COLOC="out/pQTL_vs_eQTL/pQTLs_directly_affected_colocal.csv"
if [ ! -e ${JAKOBSON_PQTL_EQTL_GENE_TARGETED_COLOC} ]; then
    echo "Extract pQTLs (Jakobson et al.) that directly affect eQTL genes (colocalized)..."
    python src/extract_directly_affected_eQTLs.py \
        --data out/pQTL_vs_eQTL/pQTLs_vs_eQTLs_colocal.csv \
        --output ${JAKOBSON_PQTL_EQTL_GENE_TARGETED_COLOC}
    echo "Done!"
fi


# Extract pQTLs that directly affect piQTL-eQTL genes
## Exact matched cases
JAKOBSON_PQTL_PIQTL_EQTL_GENE_TARGETED="out/piQTL_vs_eQTL_vs_pQTL/piQTLs_eQTLs_directly_affected_exact.csv"
if [ ! -e ${JAKOBSON_PQTL_PIQTL_EQTL_GENE_TARGETED} ]; then
    echo "Extract pQTLs (Jakobson et al.) that directly affect piQTL-eQTL genes..."
    python src/extract_directly_affected_piQTL_pQTL_eQTL.py \
        --data ${JAKOBSON_PQTL_PIQTL_EQTL_OVERLAP} \
        --output ${JAKOBSON_PQTL_PIQTL_EQTL_GENE_TARGETED}
    echo "Done!"
fi

## Co-localized cases
JAKOBSON_PQTL_PIQTL_EQTL_GENE_TARGETED_COLOC="out/piQTL_vs_eQTL_vs_pQTL/piQTLs_eQTLs_directly_affected_colocal.csv"
if [ ! -e ${JAKOBSON_PQTL_PIQTL_EQTL_GENE_TARGETED_COLOC} ]; then
    echo "Extract pQTLs (Jakobson et al.) that directly affect piQTL-eQTL genes (colocalized)..."
    python src/extract_directly_affected_piQTL_pQTL_eQTL.py \
        --data out/piQTL_vs_eQTL_vs_pQTL/piQTLs_vs_eQTLs_vs_pQTLs_colocal.csv \
        --output ${JAKOBSON_PQTL_PIQTL_EQTL_GENE_TARGETED_COLOC}
    echo "Done!"
fi
