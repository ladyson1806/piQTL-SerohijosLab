#!/usr/bin/env python3

"""
-----------------
Maps Illumina dual-indexing barcodes (N-series and S-series) to sample plate designs
and exports reference barcode files for demultiplexing pipelines.

Input:
    ../../data/pipeline/sample_position_and_barcodes.csv
        CSV with at least the columns:
            - 'N-index': N-series index ID (e.g., 'N701')
            - 'S-index': S-series index ID (e.g., 'S501')
            - 'ppi': Condition or sample identifier (e.g., 'GeneA:GeneB')

Output:
    ../../data/pipeline/PPI_reference_barcodes.csv
        Full plate design with resolved barcode sequences (CSV).
    ../../data/pipeline/PPI_reference_barcodes.txt
        Tab-separated file with columns [PPI, N_barcode, S_barcode], no header.
        Used as direct input for demultiplexing tools.

Usage:
    python 00c_barcode_reference_library.py

Dependencies:
    pandas
"""

import pandas as pd 



#### Nseries - Left inner barcode (already reverse complement)
N_series = {
    'N701': 'TCGCCTTA', 
    'N702':	'CTAGTACG',
    'N703':	'TTCTGCCT',
    'N704':	'GCTCAGGA',
    'N705':	'AGGAGTCC',
    'N706':	'CATGCCTA',
    'N707':	'GTAGAGAG',
    'N708':	'CCTCTCTG',
    'N709':	'AGCGTAGC',
    'N710':	'CAGCCTCG',
    'N711':	'TGCCTCTT',
    'N712':	'TCCTCTAC'
}

#### Sseries - Right inner barcode (already reverse complement)
S_series = {
    'S501':	'TAGATCGC',
    'S502':	'CTCTCTAT',
    'S503':	'TATCCTCT',
    'S504':	'AGAGTAGA',
    'S505':	'GTAAGGAG',
    'S506':	'ACTGCATA',
    'S507':	'AAGGAGTA',
    'S508':	'CTAAGCCT',
    'S510':	'CGTCTAAT',
    'S511':	'TCTCTCCG',
    'S513':	'TCGACTAG',
    'S515':	'TTCTAGCT'

}

def get_reference_barcode(x, series_type):
    if 'N' in series_type :
        return N_series[x]
    if 'S' in series_type :
        return S_series[x]

if __name__ == "__main__":
    plate_design = pd.read_csv('../../data/pipeline/sample_position_and_barcodes.csv')
    plate_design['N_barcode'] = plate_design['N-index'].apply(get_reference_barcode, args=('N',))
    plate_design['S_barcode'] = plate_design['S-index'].apply(get_reference_barcode, args=('S',))
    plate_design['PPI'] = [ plate_design['ppi'][i].replace(':','_') for i in plate_design.index ]
    print(plate_design)

    plate_design.to_csv('../../data/pipeline/PPI_reference_barcodes.csv', index=False)
    plate_design[['PPI', 'N_barcode', 'S_barcode']].to_csv('../../data/pipeline/PPI_reference_barcodes.txt', index=False, sep='\t', header=False)
