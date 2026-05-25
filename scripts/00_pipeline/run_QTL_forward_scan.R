library(MatrixEQTL)
library(stringr)
library(readr)
library(dplyr)
library(tibble)


args <- commandArgs(trailingOnly = TRUE)
cat(args, sep = "\n")

block_SNP_start = as.integer(args[1])
block_SNP_end = as.integer(args[2])
top_SNP = as.integer(args[3])
condition = args[4]

#### Initial Config
setwd("/home/savvy/PROJECTS/PHD/piQTL/revisions/code/")
base.dir = getwd() 

#### Genotype matrix
snps_matrix = '/../../data/genotype_information/piQTL_genotype_matrix_dec2022.txt'

#### PPI matrix
ppi_matrix = paste0('/../../results/03_ppi_estimation/logratio/all_PPI_logratio_fitness_minus_ref_delta_MTX_X_noMTX_X_for_eQTL_matrix.csv')
phe_matrix = ppi_matrix

folder = '/QTL_scan/results_0.50_final_piQTLs/'
meh_output = paste0(base.dir, folder, 'piQTL', '_', condition, '_', top_SNP, '_', 'pval_results.txt')

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


## 1. Load phenotype + genotype data into matrix
phenotype_subset = phenotype[str_detect(phenotype$Condition, condition),]
phenotype_subset <- phenotype_subset %>% column_to_rownames(var = "Condition")
phe = SlicedData$new();
phe$CreateFromMatrix(as.matrix(phenotype_subset));

snps_subset_without_hits = genotype[(genotype$snp_id >= block_SNP_start) & (genotype$snp_id <= block_SNP_end) & (genotype$snp_id != top_SNP),]
snps_subset_without_hits <- snps_subset_without_hits %>% column_to_rownames(var = "snp_id")
snps_subset_without_hits[snps_subset_without_hits == 0] <- NA
snps_block_without_hits = SlicedData$new();
snps_block_without_hits$CreateFromMatrix(as.matrix(snps_subset_without_hits));

print(snps_block_without_hits)

hits <- genotype[(genotype$snp_id == top_SNP),]
hits[hits == 0] <- NA
hits <- hits %>% column_to_rownames(var = "snp_id")
covs = SlicedData$new();
covs$CreateFromMatrix(as.matrix(hits));

print(covs)

meh_fw_scan_with_covs = Matrix_eQTL_engine(
  snps = snps_block_without_hits,
  gene = phe,
  cvrt = covs,
  output_file_name = meh_output,
  pvOutputThreshold = pvOutputThreshold,
  useModel = model,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

results_with_covariates <- meh_fw_scan_with_covs$all$eqtls
results_with_covariates$log10P <- - log10(results_with_covariates$pvalue)

## Results
## Plot the histogram of all p-values
# plot(meh);
# message('Analysis done in: ', meh$time.in.sec, ' seconds');
# message('Detected QTLs:');
# show(meh$all$eqtls);