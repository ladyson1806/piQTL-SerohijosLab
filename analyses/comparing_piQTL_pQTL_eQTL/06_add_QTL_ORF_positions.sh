#!/bin/sh

cd $(dirname ${0})

# Make a symbolic link to the SGD ORF position reference file
SGD_ORF_POSITIONS_REF="data/SGD_ORF_positions_reference.csv"
if [ ! -e ${SGD_ORF_POSITIONS_REF} ]; then
    ln -s ../../../references/SGD_database/orf_genomic_1000_R64-3-1.csv ${SGD_ORF_POSITIONS_REF}
fi


##############################################
# Add ORF positons to piQTL_vs_pQTL results
##############################################

## Exact overlaps
PIQTL_PQTL_OVERLAP_EXACT_WITH_ORF_POSITIONS="out/piQTL_vs_pQTL/piQTLs_vs_pQTLs_exact_with_ORF_positions.csv"
if [ ! -e ${PIQTL_PQTL_OVERLAP_EXACT_WITH_ORF_POSITIONS} ]; then
    echo "Add ORF positions to piQTL vs pQTL exact overlap results..."
    PIQTL_PQTL_OVERLAP_EXACT="out/piQTL_vs_pQTL/piQTLs_vs_pQTLs_exact.csv"

    python src/add_ORF_positions_to_QTL_overlap_results.py \
        --overlap_input ${PIQTL_PQTL_OVERLAP_EXACT} \
        --orf_positions_ref ${SGD_ORF_POSITIONS_REF} \
        --target_qtl_types piQTL pQTL \
        --output ${PIQTL_PQTL_OVERLAP_EXACT_WITH_ORF_POSITIONS}
    echo "Done!"
fi

## Colocalized overlaps
PIQTL_PQTL_OVERLAP_COLOCAL_WITH_ORF_POSITIONS="out/piQTL_vs_pQTL/piQTLs_vs_pQTLs_colocal_with_ORF_positions.csv"
if [ ! -e ${PIQTL_PQTL_OVERLAP_COLOCAL_WITH_ORF_POSITIONS} ]; then
    echo "Add ORF positions to piQTL vs pQTL colocalized overlap results..."
    PIQTL_PQTL_OVERLAP_COLOCAL="out/piQTL_vs_pQTL/piQTLs_vs_pQTLs_colocal.csv"
    python src/add_ORF_positions_to_QTL_overlap_results.py \
        --overlap_input ${PIQTL_PQTL_OVERLAP_COLOCAL} \
        --orf_positions_ref ${SGD_ORF_POSITIONS_REF} \
        --target_qtl_types piQTL pQTL \
        --output ${PIQTL_PQTL_OVERLAP_COLOCAL_WITH_ORF_POSITIONS}
    echo "Done!"
fi

## Exact overlaps, extracting directly affected piQTL-pQTL pairs
PIQTL_PQTL_DIRECTLY_AFFECTED_EXACT_WITH_ORF_POSITIONS="out/piQTL_vs_pQTL/piQTLs_directly_affected_exact_with_ORF_positions.csv"
if [ ! -e ${PIQTL_PQTL_DIRECTLY_AFFECTED_EXACT_WITH_ORF_POSITIONS} ]; then
    echo "Extract directly affected piQTL-pQTL pairs with ORF positions..."
    PIQTL_PQTL_DIRECTLY_AFFECTED_EXACT="out/piQTL_vs_pQTL/piQTLs_directly_affected_exact.csv"
    python src/add_ORF_positions_to_QTL_overlap_results.py \
        --overlap_input ${PIQTL_PQTL_DIRECTLY_AFFECTED_EXACT} \
        --orf_positions_ref ${SGD_ORF_POSITIONS_REF} \
        --target_qtl_types piQTL pQTL \
        --output ${PIQTL_PQTL_DIRECTLY_AFFECTED_EXACT_WITH_ORF_POSITIONS}
    echo "Done!"
fi

## Colocalized overlaps, extracting directly affected piQTL-pQTL pairs
PIQTL_PQTL_DIRECTLY_AFFECTED_COLOCAL_WITH_ORF_POSITIONS="out/piQTL_vs_pQTL/piQTLs_directly_affected_colocal_with_ORF_positions.csv"
if [ ! -e ${PIQTL_PQTL_DIRECTLY_AFFECTED_COLOCAL_WITH_ORF_POSITIONS} ]; then
    echo "Extract directly affected piQTL-pQTL pairs (colocalized) with ORF positions..."
    PIQTL_PQTL_DIRECTLY_AFFECTED_COLOCAL="out/piQTL_vs_pQTL/piQTLs_directly_affected_colocal.csv"
    python src/add_ORF_positions_to_QTL_overlap_results.py \
        --overlap_input ${PIQTL_PQTL_DIRECTLY_AFFECTED_COLOCAL} \
        --orf_positions_ref ${SGD_ORF_POSITIONS_REF} \
        --target_qtl_types piQTL pQTL \
        --output ${PIQTL_PQTL_DIRECTLY_AFFECTED_COLOCAL_WITH_ORF_POSITIONS}
    echo "Done!"
fi


##############################################
# Add ORF positons to piQTL_vs_eQTL results
##############################################

## Exact overlaps
PIQTL_EQTL_OVERLAP_EXACT_WITH_ORF_POSITIONS="out/piQTL_vs_eQTL/piQTLs_vs_eQTLs_exact_with_ORF_positions.csv"
if [ ! -e ${PIQTL_EQTL_OVERLAP_EXACT_WITH_ORF_POSITIONS} ]; then
    echo "Add ORF positions to piQTL vs eQTL exact overlap results..."
    PIQTL_EQTL_OVERLAP_EXACT="out/piQTL_vs_eQTL/piQTLs_vs_eQTLs_exact.csv"
    python src/add_ORF_positions_to_QTL_overlap_results.py \
        --overlap_input ${PIQTL_EQTL_OVERLAP_EXACT} \
        --orf_positions_ref ${SGD_ORF_POSITIONS_REF} \
        --target_qtl_types piQTL eQTL \
        --output ${PIQTL_EQTL_OVERLAP_EXACT_WITH_ORF_POSITIONS}
    echo "Done!"
fi

## Colocalized overlaps
PIQTL_EQTL_OVERLAP_COLOCAL_WITH_ORF_POSITIONS="out/piQTL_vs_eQTL/piQTLs_vs_eQTLs_colocal_with_ORF_positions.csv"
if [ ! -e ${PIQTL_EQTL_OVERLAP_COLOCAL_WITH_ORF_POSITIONS} ]; then
    echo "Add ORF positions to piQTL vs eQTL colocalized overlap results..."
    PIQTL_EQTL_OVERLAP_COLOCAL="out/piQTL_vs_eQTL/piQTLs_vs_eQTLs_colocal.csv"
    python src/add_ORF_positions_to_QTL_overlap_results.py \
        --overlap_input ${PIQTL_EQTL_OVERLAP_COLOCAL} \
        --orf_positions_ref ${SGD_ORF_POSITIONS_REF} \
        --target_qtl_types piQTL eQTL \
        --output ${PIQTL_EQTL_OVERLAP_COLOCAL_WITH_ORF_POSITIONS}
    echo "Done!"
fi

## Exact overlaps, extracting directly affected piQTL-eQTL pairs
PIQTL_EQTL_DIRECTLY_AFFECTED_EXACT_WITH_ORF_POSITIONS="out/piQTL_vs_eQTL/piQTLs_directly_affected_exact_with_ORF_positions.csv"
if [ ! -e ${PIQTL_EQTL_DIRECTLY_AFFECTED_EXACT_WITH_ORF_POSITIONS} ]; then
    echo "Extract directly affected piQTL-eQTL pairs with ORF positions..."
    PIQTL_EQTL_DIRECTLY_AFFECTED_EXACT="out/piQTL_vs_eQTL/piQTLs_directly_affected_exact.csv"
    python src/add_ORF_positions_to_QTL_overlap_results.py \
        --overlap_input ${PIQTL_EQTL_DIRECTLY_AFFECTED_EXACT} \
        --orf_positions_ref ${SGD_ORF_POSITIONS_REF} \
        --target_qtl_types piQTL eQTL \
        --output ${PIQTL_EQTL_DIRECTLY_AFFECTED_EXACT_WITH_ORF_POSITIONS}
    echo "Done!"
fi

## Colocalized overlaps, extracting directly affected piQTL-eQTL pairs
PIQTL_EQTL_DIRECTLY_AFFECTED_COLOCAL_WITH_ORF_POSITIONS="out/piQTL_vs_eQTL/piQTLs_directly_affected_colocal_with_ORF_positions.csv"
if [ ! -e ${PIQTL_EQTL_DIRECTLY_AFFECTED_COLOCAL_WITH_ORF_POSITIONS} ]; then
    echo "Extract directly affected piQTL-eQTL pairs (colocalized) with ORF positions..."
    PIQTL_EQTL_DIRECTLY_AFFECTED_COLOCAL="out/piQTL_vs_eQTL/piQTLs_directly_affected_colocal.csv"
    python src/add_ORF_positions_to_QTL_overlap_results.py \
        --overlap_input ${PIQTL_EQTL_DIRECTLY_AFFECTED_COLOCAL} \
        --orf_positions_ref ${SGD_ORF_POSITIONS_REF} \
        --target_qtl_types piQTL eQTL \
        --output ${PIQTL_EQTL_DIRECTLY_AFFECTED_COLOCAL_WITH_ORF_POSITIONS}
    echo "Done!"
fi


##############################################
# Add ORF positons to pQTL_vs_eQTL results
##############################################

## Exact overlaps
PQTL_EQTL_OVERLAP_EXACT_WITH_ORF_POSITIONS="out/pQTL_vs_eQTL/pQTLs_vs_eQTLs_exact_with_ORF_positions.csv"
if [ ! -e ${PQTL_EQTL_OVERLAP_EXACT_WITH_ORF_POSITIONS} ]; then
    echo "Add ORF positions to pQTL vs eQTL exact overlap results..."
    PQTL_EQTL_OVERLAP_EXACT="out/pQTL_vs_eQTL/pQTLs_vs_eQTLs_exact.csv"
    python src/add_ORF_positions_to_QTL_overlap_results.py \
        --overlap_input ${PQTL_EQTL_OVERLAP_EXACT} \
        --orf_positions_ref ${SGD_ORF_POSITIONS_REF} \
        --target_qtl_types pQTL eQTL \
        --output ${PQTL_EQTL_OVERLAP_EXACT_WITH_ORF_POSITIONS}
    echo "Done!"
fi

## Colocalized overlaps
PQTL_EQTL_OVERLAP_COLOCAL_WITH_ORF_POSITIONS="out/pQTL_vs_eQTL/pQTLs_vs_eQTLs_colocal_with_ORF_positions.csv"
if [ ! -e ${PQTL_EQTL_OVERLAP_COLOCAL_WITH_ORF_POSITIONS} ]; then
    echo "Add ORF positions to pQTL vs eQTL colocalized overlap results..."
    PQTL_EQTL_OVERLAP_COLOCAL="out/pQTL_vs_eQTL/pQTLs_vs_eQTLs_colocal.csv"
    python src/add_ORF_positions_to_QTL_overlap_results.py \
        --overlap_input ${PQTL_EQTL_OVERLAP_COLOCAL} \
        --orf_positions_ref ${SGD_ORF_POSITIONS_REF} \
        --target_qtl_types pQTL eQTL \
        --output ${PQTL_EQTL_OVERLAP_COLOCAL_WITH_ORF_POSITIONS}
    echo "Done!"
fi

## Exact overlaps, extracting directly affected pQTL-eQTL pairs
PQTL_EQTL_DIRECTLY_AFFECTED_EXACT_WITH_ORF_POSITIONS="out/pQTL_vs_eQTL/pQTLs_directly_affected_exact_with_ORF_positions.csv"
if [ ! -e ${PQTL_EQTL_DIRECTLY_AFFECTED_EXACT_WITH_ORF_POSITIONS} ]; then
    echo "Extract directly affected pQTL-eQTL pairs with ORF positions..."
    PQTL_EQTL_DIRECTLY_AFFECTED_EXACT="out/pQTL_vs_eQTL/pQTLs_directly_affected_exact.csv"
    python src/add_ORF_positions_to_QTL_overlap_results.py \
        --overlap_input ${PQTL_EQTL_DIRECTLY_AFFECTED_EXACT} \
        --orf_positions_ref ${SGD_ORF_POSITIONS_REF} \
        --target_qtl_types pQTL eQTL \
        --output ${PQTL_EQTL_DIRECTLY_AFFECTED_EXACT_WITH_ORF_POSITIONS}
    echo "Done!"
fi

## Colocalized overlaps, extracting directly affected pQTL-eQTL pairs
PQTL_EQTL_DIRECTLY_AFFECTED_COLOCAL_WITH_ORF_POSITIONS="out/pQTL_vs_eQTL/pQTLs_directly_affected_colocal_with_ORF_positions.csv"
if [ ! -e ${PQTL_EQTL_DIRECTLY_AFFECTED_COLOCAL_WITH_ORF_POSITIONS} ]; then
    echo "Extract directly affected pQTL-eQTL pairs (colocalized) with ORF positions..."
    PQTL_EQTL_DIRECTLY_AFFECTED_COLOCAL="out/pQTL_vs_eQTL/pQTLs_directly_affected_colocal.csv"
    python src/add_ORF_positions_to_QTL_overlap_results.py \
        --overlap_input ${PQTL_EQTL_DIRECTLY_AFFECTED_COLOCAL} \
        --orf_positions_ref ${SGD_ORF_POSITIONS_REF} \
        --target_qtl_types pQTL eQTL \
        --output ${PQTL_EQTL_DIRECTLY_AFFECTED_COLOCAL_WITH_ORF_POSITIONS}
    echo "Done!"
fi


##############################################
# Add ORF positons to piQTL_vs_eQTL_vs_pQTL results
##############################################

## Exact overlaps
PIQTL_EQTL_PQTL_OVERLAP_EXACT_WITH_ORF_POSITIONS="out/piQTL_vs_eQTL_vs_pQTL/piQTLs_vs_eQTLs_vs_pQTLs_exact_with_ORF_positions.csv"
if [ ! -e ${PIQTL_EQTL_PQTL_OVERLAP_EXACT_WITH_ORF_POSITIONS} ]; then
    echo "Add ORF positions to piQTL vs eQTL vs pQTL exact overlap results..."
    PIQTL_EQTL_PQTL_OVERLAP_EXACT="out/piQTL_vs_eQTL_vs_pQTL/piQTLs_vs_eQTLs_vs_pQTLs_exact.csv"
    python src/add_ORF_positions_to_QTL_overlap_results.py \
        --overlap_input ${PIQTL_EQTL_PQTL_OVERLAP_EXACT} \
        --orf_positions_ref ${SGD_ORF_POSITIONS_REF} \
        --target_qtl_types piQTL eQTL pQTL \
        --output ${PIQTL_EQTL_PQTL_OVERLAP_EXACT_WITH_ORF_POSITIONS}
    echo "Done!"
fi

## Colocalized overlaps
PIQTL_EQTL_PQTL_OVERLAP_COLOCAL_WITH_ORF_POSITIONS="out/piQTL_vs_eQTL_vs_pQTL/piQTLs_vs_eQTLs_vs_pQTLs_colocal_with_ORF_positions.csv"
if [ ! -e ${PIQTL_EQTL_PQTL_OVERLAP_COLOCAL_WITH_ORF_POSITIONS} ]; then
    echo "Add ORF positions to piQTL vs eQTL vs pQTL colocalized overlap results..."
    PIQTL_EQTL_PQTL_OVERLAP_COLOCAL="out/piQTL_vs_eQTL_vs_pQTL/piQTLs_vs_eQTLs_vs_pQTLs_colocal.csv"
    python src/add_ORF_positions_to_QTL_overlap_results.py \
        --overlap_input ${PIQTL_EQTL_PQTL_OVERLAP_COLOCAL} \
        --orf_positions_ref ${SGD_ORF_POSITIONS_REF} \
        --target_qtl_types piQTL eQTL pQTL \
        --output ${PIQTL_EQTL_PQTL_OVERLAP_COLOCAL_WITH_ORF_POSITIONS}
    echo "Done!"
fi

## Extract directly affected piQTL-eQTL-pQTL pairs with ORF positions
PIQTL_EQTL_PQTL_DIRECTLY_AFFECTED_EXACT_WITH_ORF_POSITIONS="out/piQTL_vs_eQTL_vs_pQTL/piQTLs_eQTLs_directly_affected_exact_with_ORF_positions.csv"
if [ ! -e ${PIQTL_EQTL_PQTL_DIRECTLY_AFFECTED_EXACT_WITH_ORF_POSITIONS} ]; then
    echo "Extract directly affected piQTL-eQTL-pQTL pairs with ORF positions..."
    PIQTL_EQTL_PQTL_DIRECTLY_AFFECTED_EXACT="out/piQTL_vs_eQTL_vs_pQTL/piQTLs_eQTLs_directly_affected_exact.csv"
    python src/add_ORF_positions_to_QTL_overlap_results.py \
        --overlap_input ${PIQTL_EQTL_PQTL_DIRECTLY_AFFECTED_EXACT} \
        --orf_positions_ref ${SGD_ORF_POSITIONS_REF} \
        --target_qtl_types piQTL eQTL pQTL \
        --output ${PIQTL_EQTL_PQTL_DIRECTLY_AFFECTED_EXACT_WITH_ORF_POSITIONS}
    echo "Done!"
fi

## Extract directly affected piQTL-eQTL-pQTL pairs (colocalized) with ORF positions
PIQTL_EQTL_PQTL_DIRECTLY_AFFECTED_COLOCAL_WITH_ORF_POSITIONS="out/piQTL_vs_eQTL_vs_pQTL/piQTLs_eQTLs_directly_affected_colocal_with_ORF_positions.csv"
if [ ! -e ${PIQTL_EQTL_PQTL_DIRECTLY_AFFECTED_COLOCAL_WITH_ORF_POSITIONS} ]; then
    echo "Extract directly affected piQTL-eQTL-pQTL pairs (colocalized) with ORF positions..."
    PIQTL_EQTL_PQTL_DIRECTLY_AFFECTED_COLOCAL="out/piQTL_vs_eQTL_vs_pQTL/piQTLs_eQTLs_directly_affected_colocal.csv"
    python src/add_ORF_positions_to_QTL_overlap_results.py \
        --overlap_input ${PIQTL_EQTL_PQTL_DIRECTLY_AFFECTED_COLOCAL} \
        --orf_positions_ref ${SGD_ORF_POSITIONS_REF} \
        --target_qtl_types piQTL eQTL pQTL \
        --output ${PIQTL_EQTL_PQTL_DIRECTLY_AFFECTED_COLOCAL_WITH_ORF_POSITIONS}
    echo "Done!"
fi
