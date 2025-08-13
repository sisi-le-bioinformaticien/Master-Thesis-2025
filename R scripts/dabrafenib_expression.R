# --- Libraries ---
library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(broom)
library(purrr)
library(readr)
library(viridis)


# Load data
rna_file <- "./rnaseq_merged_20250117/rnaseq_merged_20250117_filtered.csv"
rna_data <- read.csv(rna_file, sep = ",") %>%
  mutate(htseq_read_count = as.integer(htseq_read_count))

gdsc1 <- read_excel("./data/GDSC1_fitted_dose_response_27Oct23.xlsx")
gdsc2 <- read_excel("./data/GDSC2_fitted_dose_response_27Oct23.xlsx")
gdsc <- bind_rows(gdsc1, gdsc2)

braf_mut_file <- "./BRAF_cell_line_mut_index2.csv"
braf_status <- read.csv(braf_mut_file) %>%
  select(model_name, BRAF_mutations) %>%
  mutate(condition = case_when(
    grepl("V600E", BRAF_mutations, ignore.case = TRUE) ~ "BRAF_V600E",
    grepl("Wildtype", BRAF_mutations, ignore.case = TRUE) ~ "BRAF_WT",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(condition))
colnames(braf_status)[colnames(braf_status) == "model_name"] <- "CELL_LINE_NAME"

# Filter RNA-seq data to relevant cell lines 
selected_lines <- unique(braf_status$CELL_LINE_NAME)
expr_data <- rna_data %>%
  filter(model_name %in% selected_lines) %>%
  group_by(model_name, gene_symbol) %>%
  summarise(htseq_read_count = sum(htseq_read_count), .groups = "drop")

# Build expression matrix 
expr_matrix <- expr_data %>%
  pivot_wider(names_from = model_name, values_from = htseq_read_count) %>%
  column_to_rownames("gene_symbol") %>%
  as.matrix()

expr_matrix <- expr_matrix[rowSums(expr_matrix) > 10, ]  # Filter low-count genes

# Create metadata for mutation status
meta <- braf_status %>%
  filter(CELL_LINE_NAME %in% colnames(expr_matrix)) %>%
  distinct(CELL_LINE_NAME, condition) %>%
  column_to_rownames("CELL_LINE_NAME")

# Align expression and metadata
expr_matrix <- expr_matrix[, rownames(meta)]

# DESeq2 Analysis 
dds <- DESeqDataSetFromMatrix(countData = expr_matrix, colData = meta, design = ~condition)
dds <- DESeq(dds)
res <- results(dds)

# Top differentially expressed genes
res <- as.data.frame(res)
res$gene <- rownames(res)

top_genes2 <- res %>%
  filter(!is.na(padj)) %>%
  arrange(-abs(log2FoldChange)) %>%
  slice_head(n = 12) 

print(top_genes2)



# PART 2: Compute correlation to GDSC data with boltz predictions


df <- read_csv("drug_comparison_IC50.csv")

colnames(vst_mat) <- toupper(gsub("-", "", colnames(vst_mat)))
gdsc$CELL_LINE_NAME <- toupper(gsub("-", "", gdsc$CELL_LINE_NAME))

gdsc$DRUG_NAME <- gsub("[-\\s()]", "", tolower(gdsc$DRUG_NAME))

print(rownames(top_genes2))
top_genes <- intersect(rownames(top_genes2), rownames(vst_mat))


expr_score_for_drug <- function(drug, genes, vst_mat, gdsc) {
  tested_lines <- gdsc %>%
    filter(DRUG_NAME == drug) %>%
    pull(CELL_LINE_NAME) %>%
    unique()
  
  matching_lines <- intersect(colnames(vst_mat), tested_lines)
  if (length(matching_lines) == 0) return(NA_real_)
  
  expr_subset <- vst_mat[genes, matching_lines, drop = FALSE]
  median(expr_subset, na.rm = TRUE)
}

df$expr_score_perdrug <- map_dbl(df$Drug, expr_score_for_drug,
                                 genes = top_genes,
                                 vst_mat = vst_mat,
                                 gdsc = gdsc)

df2 <- df %>% filter(!is.na(expr_score_perdrug))

m2 <- lm(GDSC_IC50 ~ Boltz_IC50_BRAF + expr_score_perdrug , data = df2)


# Spearman correlation 
spearman_m2 <- cor(predict(m2), df2$GDSC_IC50, method = "spearman")

cat("Spearman (Model 2 vs GDSC):", round(spearman_m2, 3), "\n")

# Plot predicted vs observed
df2$predicted_m2 <- predict(m1)

ggplot(df2, aes(x = predicted_m2, y = GDSC_IC50)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed") +
  labs(
    title = "Correlation between GDSC and Boltz2 + Expression Score",
    subtitle = paste0("Spearman =", round(spearman_m2, 3)),
    x = "Predicted log(IC50) (combined model)",
    y = "Observed GDSC log(IC50)"
  ) +
  theme_minimal(base_size = 13)

write_csv(df2, "drug_comparison_IC50_with_expr.csv")




