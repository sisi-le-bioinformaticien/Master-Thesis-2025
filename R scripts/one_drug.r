


library(readxl)
library(ggplot2)
library(patchwork)
library(tidyr)
library(dplyr)
library(DESeq2)
library(pheatmap)



SMALL_DATA <- read.csv("./rnaseq_merged_20250117/rnaseq_merged_20250117_filtered.csv",sep=",") 


# DRUG Cisplatin|Fluorouracil|Rapamycin|Methotrexate|Flavopiridol|Irinotecan


data1 <- read_excel("./data/GDSC1_fitted_dose_response_27Oct23.xlsx") 
data2 <- read_excel("./data/GDSC2_fitted_dose_response_27Oct23.xlsx") 
combined_data <- bind_rows(data1, data2)

data_cisplatin <- combined_data %>%
  filter(grepl("GNF-2", DRUG_NAME)) %>%
  select(CELL_LINE_NAME,Z_SCORE)

resistant_cell_lines <- data_cisplatin %>%
  filter(Z_SCORE > 2) %>%
  select(CELL_LINE_NAME) %>%
  pull(CELL_LINE_NAME)


sensible_cell_lines <- data_cisplatin %>%
  filter(Z_SCORE < -2) %>%
  select(CELL_LINE_NAME) %>%
  pull(CELL_LINE_NAME)


DRUG_SPECIFIC_DATA <- SMALL_DATA %>%
  filter(model_name %in% sensible_cell_lines | model_name %in% resistant_cell_lines )




EXPRESSION_MATRIX <- DRUG_SPECIFIC_DATA %>%
  select(model_name, gene_symbol, htseq_read_count) %>%
  pivot_wider(
    names_from = model_name,
    values_from = htseq_read_count
    
  )


EXPRESSION_MATRIX <- as.data.frame(EXPRESSION_MATRIX)
rownames(EXPRESSION_MATRIX) <- EXPRESSION_MATRIX[[1]]  
EXPRESSION_MATRIX <- EXPRESSION_MATRIX[ , -1] 
EXPRESSION_MATRIX[] <- lapply(EXPRESSION_MATRIX, function(x) unlist(x))
EXPRESSION_MATRIX <- as.matrix(EXPRESSION_MATRIX)
dim(EXPRESSION_MATRIX)


condition <- ifelse(
  colnames(EXPRESSION_MATRIX) %in% resistant_cell_lines, "resistant",
  ifelse(colnames(EXPRESSION_MATRIX) %in% sensible_cell_lines, "sensitive", NA)
)
condition <- factor(condition, levels = c("resistant", "sensitive"))

metadata <- data.frame(condition)
rownames(metadata) <- colnames(EXPRESSION_MATRIX)


dds <- DESeqDataSetFromMatrix(
  countData = EXPRESSION_MATRIX,
  colData = metadata,
  design = ~ condition
)

dds <- DESeq(dds)
res <- results(dds)


# PCA plot
vsd <- vst(dds, blind = TRUE)  

# Heatmap 
heatmap_data <- assay(vsd)

row_vars <- apply(heatmap_data, 1, var)
sum(row_vars == 0) 
heatmap_data_filtered <- heatmap_data[row_vars > 0, ]

condition_colors <- list(condition = c(resistant = "orange", sensitive = "green"))

pheatmap(heatmap_data_filtered,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         show_rownames = FALSE,
         show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         annotation_col = metadata,
         annotation_colors = condition_colors
)
print("done")

top_genes <- head(res[order(res$padj), ], 20)

heatmap_top_genes <- heatmap_data_filtered[rownames(top_genes), ]

pheatmap(heatmap_top_genes,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         annotation_col = metadata,
         annotation_colors = condition_colors
)


