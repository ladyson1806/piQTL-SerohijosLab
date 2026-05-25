library(MatrixEQTL)
library(stringr)
library(readr)
library(dplyr)
library(tibble)

# Running this script on the cluster with:
# export MKL_THREADING_LAYER=sequential 

#### Initial Config

base.dir = getwd() 

#### Genotype matrix
snps_matrix = '/../../data/genotype_information/piQTL_genotype_matrix_dec2022.txt'

#### PPI matrix
# We have several versions of the PPI matrix, depending on how we want to define the PPI changes (before downsampling, normalized or not). We can choose which one to use for the eQTL mapping.
# ppi_matrix = paste0('/../../results/03_ppi_estimation/logratio/all_PPI_logratio_fitness_before_downsampling_minus_ref_for_eQTL_matrix.csv')
# ppi_matrix = paste0('/../../results/03_ppi_estimation/logratio/all_PPI_logratio_fitness_minus_ref_delta_MTX_X_noMTX_X_for_eQTL_matrix.csv')
ppi_matrix = paste0('/../../results/03_ppi_estimation/logratio/all_PPI_logratio_fitness_minus_ref_delta_MTX_X_noMTX_X_normalized_for_eQTL_matrix.csv')
phe_matrix = ppi_matrix

folder = '/../../results/04_eQTL_matrix/'
meh_output = paste0(base.dir, folder, 'piQTL', '_', 'june2024', '_', 'pval_results.txt')
print(paste0('Dump QTL mapping results to: ', meh_output))

#### Parameters

SNP_file_name = paste0(base.dir, snps_matrix)
genotype = readr::read_csv(SNP_file_name)
phe_file_name = paste0(base.dir, phe_matrix)
phenotype = readr::read_csv(phe_file_name)
## Define the model
model = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
## Threshold for p-values
pvOutputThreshold = 1;
# Error covariance matrix
errorCovariance = numeric(); # Set to numeric() for identity.

## 1. Load genotype data into matrix
snps = SlicedData$new();
snps$fileDelimiter = ",";
snps$fileOmitCharacters = "0"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## 2. Load phenotypic measurement data into matrix
phe = SlicedData$new();
phe$fileDelimiter = ",";
phe$fileOmitCharacters = "NA"; # denote missing values;
phe$fileSkipRows = 1;          # one row of column labels
phe$fileSkipColumns = 1;       # one column of row labels
phe$fileSliceSize = 2000;      # read file in slices of 2,000 rows
phe$LoadFile(phe_file_name);

## 3. Run the analysis
meh = Matrix_eQTL_engine(
snps = snps,
  gene = phe,
  output_file_name = meh_output,
  pvOutputThreshold = pvOutputThreshold,
  useModel = model,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

## Results
message('Analysis done in: ', meh$time.in.sec, ' seconds');