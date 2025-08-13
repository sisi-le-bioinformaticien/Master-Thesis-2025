

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install(c("GSVA", "GSEABase"))
#BiocManager::install("progeny")

#install.packages("https://cran.r-project.org/src/contrib/Archive/matrixStats/matrixStats_0.63.0.tar.gz", 
#                 repos = NULL, 
#                 type = "source")

#packageVersion("matrixStats")


#import libs
library(readxl)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(GSVA)
library(GSEABase)
library(tidyr)
library(dplyr)

#browseVignettes("GSVA")

#ALL_DATA <- read.csv("./data/rnaseq_merged_20250117.csv",sep=",") 
#SMALL_DATA <- ALL_DATA %>% 
#  select(model_name, gene_symbol,ensembl_gene_id, rsem_tpm) %>% 
#  filter(if_any(everything(), ~ !is.na(rsem_tpm)))
#write.csv(SMALL_DATA, "./data/rnaseq_merged_20250117_filtered_tpm.csv", row.names = FALSE)


#import files
SMALL_DATA_TPM <- read.csv("./data/rnaseq_merged_20250117_filtered_tpm.csv",sep=",") 
gene_sets <- getGmt("./data/c6.all.v2024.1.Hs.symbols.gmt")
data1 <- read_excel("./data/GDSC1_fitted_dose_response_27Oct23.xlsx") 
data2 <- read_excel("./data/GDSC2_fitted_dose_response_27Oct23.xlsx") 
combined_data <- bind_rows(data1, data2)


gene_count_per_model <- table(SMALL_DATA_TPM$model_name)
table(gene_count_per_model)



first_100_models <- SMALL_DATA_TPM %>%
  distinct(model_name) %>%
  slice(1:100) %>%
  pull(model_name)

EXPRESSION_MATRIX_INIT <- SMALL_DATA_TPM %>%
  filter(model_name %in% first_100_models)




EXPRESSION_MATRIX_INIT <- EXPRESSION_MATRIX_INIT %>%
  group_by(gene_symbol, model_name) %>%
  summarise(rsem_tpm = mean(rsem_tpm), .groups = "drop")%>%
  pivot_wider(
    names_from = model_name,
    values_from = rsem_tpm
  )


EXPRESSION_MATRIX <- as.data.frame(EXPRESSION_MATRIX_INIT)
rownames(EXPRESSION_MATRIX) <- EXPRESSION_MATRIX$gene_symbol
EXPRESSION_MATRIX$gene_symbol <- NULL
EXPRESSION_MATRIX <- as.matrix(EXPRESSION_MATRIX)
EXPRESSION_MATRIX <- log2(EXPRESSION_MATRIX + 1)

dim(EXPRESSION_MATRIX)




# compute gsva scores
gsva_scores <- gsva(EXPRESSION_MATRIX, gene_sets, method = "gsva", min.sz = 5, kcdf = "Gaussian")

# TEST other lib for pathway scrores
#BiocManager::install("progeny")





# Filter datasets to keep common cell lines only
cell_lines_gsva <- colnames(gsva_scores)
cell_lines_ic50 <- unique(combined_data$CELL_LINE_NAME)
common_cell_lines <- intersect(cell_lines_gsva, cell_lines_ic50)

gsva_common <- gsva_scores[, common_cell_lines]
ic50_common <- combined_data %>%
  filter(CELL_LINE_NAME %in% common_cell_lines)
ic50_avg <- ic50_common %>%
  group_by(DRUG_NAME, CELL_LINE_NAME) %>%
  summarize(LN_IC50 = mean(LN_IC50, na.rm = TRUE), .groups = "drop")


ic50_matrix <- ic50_avg %>%
  select(DRUG_NAME, CELL_LINE_NAME, LN_IC50) %>%
  pivot_wider(names_from = CELL_LINE_NAME, values_from = LN_IC50)

ic50_matrix <- as.data.frame(ic50_matrix)
rownames(ic50_matrix) <- ic50_matrix$DRUG_NAME
ic50_matrix$DRUG_NAME <- NULL

# Convert to matrix
ic50_matrix <- as.matrix(ic50_matrix)

common_cols <- intersect(colnames(gsva_common), colnames(ic50_matrix))
gsva_common_final <- gsva_common[, common_cols]
ic50_common_final <- ic50_matrix[, common_cols]

table(sapply(as.data.frame(ic50_common_final), class))



# Compute correlations
cor_matrix <- cor(t(gsva_common), t(ic50_common_final), method = "spearman", use = "pairwise.complete.obs")

sum(is.na(cor_matrix))



which(abs(cor_matrix) > 0.5, arr.ind = TRUE) 

color_scale <- colorRampPalette(c("blue", "white", "red"))(100) 
pheatmap(cor_matrix,
         cluster_rows = F,
         cluster_cols = F,
         #scale = "row",
         show_rownames = T,
         show_colnames = T,
         drop=T,
         fontsize = 4,
         color = color_scale,
         legend= T
)
print("done")

significant_corrs <- which(abs(cor_matrix) > 0.5, arr.ind = TRUE)
print(paste("Number of significant correlations (> 0.5):", nrow(significant_corrs)))

# Summary of correlation values
summary(cor_matrix)




### PROGENY

library(progeny)
progeny_scores <- progeny(EXPRESSION_MATRIX, scale=TRUE, organism="Human", top=100)
progeny_scores <-t(progeny_scores)

cell_lines_progeny <- colnames(progeny_scores)
cell_lines_ic50 <- unique(combined_data$CELL_LINE_NAME)
common_cell_lines <- intersect(cell_lines_progeny, cell_lines_ic50)

progeny_common <- progeny_scores[, common_cell_lines]
ic50_common <- combined_data %>%
  filter(CELL_LINE_NAME %in% common_cell_lines)
ic50_avg <- ic50_common %>%
  group_by(DRUG_NAME, CELL_LINE_NAME) %>%
  summarize(LN_IC50 = mean(LN_IC50, na.rm = TRUE), .groups = "drop")


ic50_matrix <- ic50_avg %>%
  select(DRUG_NAME, CELL_LINE_NAME, LN_IC50) %>%
  pivot_wider(names_from = CELL_LINE_NAME, values_from = LN_IC50)

ic50_matrix <- as.data.frame(ic50_matrix)
rownames(ic50_matrix) <- ic50_matrix$DRUG_NAME
ic50_matrix$DRUG_NAME <- NULL

# Convert to matrix
ic50_matrix <- as.matrix(ic50_matrix)

common_cols <- intersect(colnames(progeny_common), colnames(ic50_matrix))
progeny_common_final <- progeny_common[, common_cols]
ic50_common_final <- ic50_matrix[, common_cols]

table(sapply(as.data.frame(ic50_common_final), class))


# Compute correlations
cor_matrix <- cor(t(progeny_common), t(ic50_common_final), method = "spearman", use = "pairwise.complete.obs")


hist(which(abs(cor_matrix) > 0.5, arr.ind = TRUE), xlab = rownames(cor_matrix) )

color_scale <- colorRampPalette(c("blue", "white", "red"))(100) 
pheatmap(cor_matrix,
         cluster_rows = F,
         cluster_cols = F,
         #scale = "row",
         show_rownames = T,
         show_colnames = T,
         drop=T,
         fontsize = 4,
         color = color_scale,
         legend= T
)
print("done")
print("done")
