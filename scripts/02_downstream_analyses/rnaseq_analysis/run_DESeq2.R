library(DESeq2)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(data.table)
library(AnnotationDbi)
library(org.Sc.sgd.db)

file_path = getwd()

#### Inputs
sample_info = fread('rna_seq_metadata.csv')
sample_info$DRUG <- relevel(as.factor(sample_info$DRUG), "NoDrug")

featurecounts = fread('average_raw_counts.csv')

#### Data preprocessing
gene_count_matrix <- as.data.frame(featurecounts) %>% 
  column_to_rownames("Geneid") %>% # turn the geneid column into rownames
  as.matrix()

keep <- rowSums(gene_count_matrix) > 0
gene_count_matrix <- gene_count_matrix[keep,]
paste0('Number of kept genes:', nrow(gene_count_matrix))


########PCA exploration########
dds.raw <- DESeqDataSetFromMatrix(countData = round(gene_count_matrix),
                                     colData = sample_info,
                                     design = ~ Strain_ID + DRUG)

vstcounts <- vst(dds.raw)
write.table(assay(vstcounts), './normalized_counts.csv', sep=",", row.names=T)

pca_data = plotPCA(vstcounts, ntop=5000, intgroup=c("Strain_ID","DRUG"), returnData=TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(x= PC1, y = PC2)) +
  geom_point(size= 3, aes(color=DRUG)) +
  geom_text_repel(size= 4, aes(label=sample_info$Strain_ID)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  scale_color_manual(values = c("deepskyblue2", "azure4"), name = "DRUG", labels = c("Fluconazole", "No Drug")) +
  ggtitle("PCA NoDrug vs Fluconazole") +
  theme_light()

dds.diff = DESeq(dds.raw)
# compute the results object and set the alpha for the FDR at the threshold that I want to choose after
res = results(dds.diff, contrast=c("DRUG", "Fluconazole", "NoDrug"), alpha = 0.05)
all <- res %>% as.data.frame()
dim(all)
# significant genes
sign = res %>% as.data.frame() %>%  filter(padj < 0.05)
write.table(sign,"./all_differently_expressed_genes.csv", quote=F, sep=",", col.names = NA)
dim(sign)

lfc_apeglm = lfcShrink(dds.diff, res = res, coef = "DRUG_Fluconazole_vs_NoDrug" , type="apeglm")
#add information that needs to be displayed in the plot
apeglm_for_plot = lfc_apeglm %>% as.data.frame() %>%  
  rownames_to_column(var = "locus_id") %>% dplyr::select(c("locus_id", "log2FoldChange", "padj")) %>% na.omit() %>% 
  mutate(significance = ifelse(padj < 0.05, "signif", "non_signif")) %>% 
  mutate(symbol = mapIds(org.Sc.sgd.db, 
                         keys= .$locus_id,
                         column="GENENAME", 
                         keytype="ORF",
                         multiVals="first")) %>% 
  mutate(annot = ifelse(log2FoldChange < -2, symbol, ifelse(log2FoldChange > 2, symbol, ""))) %>% 
  mutate(annot = ifelse(is.na(annot), locus_id, annot))

# plot, label for the genes that have a log2FC > |2|
ggplot(apeglm_for_plot, aes(log2FoldChange, -log10(padj))) + 
  geom_point(aes(col = significance), alpha = 0.5) + theme_light() +
  scale_color_manual(name = "FDR 5%", labels = c("Non significant", "Significant"), values = c("lightgrey", "red")) +
  geom_text_repel(aes(label = annot), size = 4) + 
  ylab("-log10(FDR)") + xlab("log2FoldChange Fluconazole vs. NoDrug") +
  theme(legend.position = "bottom", legend.direction = "horizontal")
