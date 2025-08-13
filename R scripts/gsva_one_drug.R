library(readxl)
library(dplyr)
library(tidyr)
library(DESeq2)
library(pheatmap)
library(GSVA)
library(GSEABase)
library(limma)
library(ggplot2)

# Load data
SMALL_DATA <- read.csv("./data/rnaseq_merged_20250117_filtered.csv", sep=",")
SMALL_DATA_TPM <- read.csv("./data/rnaseq_merged_20250117_filtered_tpm.csv", sep=",") 
gene_sets <- getGmt("./data/c6.all.v2024.1.Hs.symbols.gmt")

# Load drug response data
data1 <- read_excel("./data/GDSC1_fitted_dose_response_27Oct23.xlsx") 
data2 <- read_excel("./data/GDSC2_fitted_dose_response_27Oct23.xlsx") 
combined_data <- bind_rows(data1, data2)



drug_name <- "Methotrexate" 



data_drug <- combined_data %>%
  filter(DRUG_NAME == drug_name) %>%
  select(CELL_LINE_NAME, Z_SCORE)

# Identify resistant and sensitive cell lines
resistant_cell_lines <- data_drug %>%
  filter(Z_SCORE > 2) %>%
  pull(CELL_LINE_NAME)

sensible_cell_lines <- data_drug %>%
  filter(Z_SCORE < -2) %>%
  pull(CELL_LINE_NAME)

cat("Found", length(resistant_cell_lines), "resistant and", 
    length(sensible_cell_lines), "sensitive cell lines for", drug_name, "\n")

if (length(resistant_cell_lines) < 3 | length(sensible_cell_lines) < 3) {
  stop("Not enough resistant or sensitive cell lines to perform analysis")
}

# Differential expression analysis
DRUG_SPECIFIC_DATA <- SMALL_DATA %>%
  filter(model_name %in% sensible_cell_lines | model_name %in% resistant_cell_lines)

EXPRESSION_MATRIX <- DRUG_SPECIFIC_DATA %>%
  select(model_name, gene_symbol, htseq_read_count) %>%
  pivot_wider(names_from = model_name, values_from = htseq_read_count) %>%
  as.data.frame()

rownames(EXPRESSION_MATRIX) <- EXPRESSION_MATRIX[[1]]  
EXPRESSION_MATRIX <- EXPRESSION_MATRIX[, -1] 
EXPRESSION_MATRIX[] <- lapply(EXPRESSION_MATRIX, function(x) unlist(x))
EXPRESSION_MATRIX <- as.matrix(EXPRESSION_MATRIX)

# Create condition vector
condition <- ifelse(
  colnames(EXPRESSION_MATRIX) %in% resistant_cell_lines, "resistant",
  ifelse(colnames(EXPRESSION_MATRIX) %in% sensible_cell_lines, "sensitive", NA)
)
condition <- factor(condition, levels = c("resistant", "sensitive"))

metadata <- data.frame(condition)
rownames(metadata) <- colnames(EXPRESSION_MATRIX)

# Run DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = EXPRESSION_MATRIX,
  colData = metadata,
  design = ~ condition
)

dds <- DESeq(dds)
res <- results(dds)

# PCA plot
vsd <- vst(dds, blind = TRUE)
pca_plot <- plotPCA(vsd, intgroup = "condition") + 
  ggtitle(paste("PCA plot for", drug_name, "response"))
print(pca_plot)

# Heatmap of top differentially expressed genes
top_genes <- head(res[order(res$padj), ], 20)
heatmap_data <- assay(vsd)
heatmap_top_genes <- heatmap_data[rownames(top_genes), ]

pheatmap(heatmap_top_genes,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         annotation_col = metadata,
         main = paste("Top DE genes for", drug_name)
)

# GSVA pathway analysis
# Prepare TPM data
TPM_MATRIX <- SMALL_DATA_TPM %>%
  filter(model_name %in% c(resistant_cell_lines, sensible_cell_lines)) %>%
  group_by(gene_symbol, model_name) %>%
  summarise(rsem_tpm = mean(rsem_tpm), .groups = "drop") %>%
  pivot_wider(names_from = model_name, values_from = rsem_tpm) %>%
  as.data.frame()

rownames(TPM_MATRIX) <- TPM_MATRIX$gene_symbol
TPM_MATRIX$gene_symbol <- NULL
TPM_MATRIX <- as.matrix(TPM_MATRIX)
TPM_MATRIX <- log2(TPM_MATRIX + 1)

# Compute GSVA scores
gsva_scores <- gsva(TPM_MATRIX, gene_sets, method = "gsva", min.sz = 5, kcdf = "Gaussian")

# Compare pathway activity between resistant and sensitive
pathway_condition <- ifelse(
  colnames(gsva_scores) %in% resistant_cell_lines, "resistant", "sensitive"
)

# Differential pathway analysis
pathway_design <- model.matrix(~ pathway_condition)
fit <- lmFit(gsva_scores, pathway_design)
fit <- eBayes(fit)
pathway_results <- topTable(fit, number = Inf, adjust.method = "BH")

cat("\nTop differentially active pathways:\n")
print(head(pathway_results[c("logFC", "adj.P.Val")], 10))

# Heatmap of significant pathways
sig_pathways <- rownames(pathway_results)[pathway_results$adj.P.Val < 0.1]
if (length(sig_pathways) > 1) {
  pathway_metadata <- data.frame(Response = pathway_condition)
  rownames(pathway_metadata) <- colnames(gsva_scores)
  
  pheatmap(gsva_scores[sig_pathways, ],
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           scale = "row",
           show_rownames = TRUE,
           show_colnames = FALSE,
           annotation_col = pathway_metadata,
           main = paste("Pathway activity for", drug_name),
           fontsize_row = 6
  )
} else {
  cat("No significantly differentially active pathways found at FDR < 0.1\n")
}
