# Load required libraries
library(dplyr)
library(readxl)
library(tidyr)
library(tibble)

# Load expression data and drug response data
MUT_DATA <- read.csv("./data/mutations_all.csv", sep=",") 
SMALL_DATA <- read.csv("./data/rnaseq_merged_20250117_filtered.csv",sep=",") 
data1 <- read_excel("./data/GDSC1_fitted_dose_response_27Oct23.xlsx") 
data2 <- read_excel("./data/GDSC2_fitted_dose_response_27Oct23.xlsx") 
combined_data <- bind_rows(data1, data2)


# 1. Prepare expression matrix (genes x cell lines)
expr_wide <- SMALL_DATA %>%
  select(model_name, ensembl_gene_id, htseq_read_count) %>%
  pivot_wider(names_from = model_name, values_from = htseq_read_count)

expr_wide[-1] <- lapply(expr_wide[-1], function(x) as.numeric(as.character(x)))

expr_mat <- expr_wide %>%
  column_to_rownames("ensembl_gene_id") %>%
  as.matrix()


expr_mat <- log2(expr_mat + 1) 

expr_mat_ranked <- t(apply(expr_mat, 1, rank, ties.method = "average"))

# 2. Prepare IC50 matrix (drugs x cell lines)
ic50_summary <- combined_data %>%
  group_by(DRUG_NAME, CELL_LINE_NAME) %>%
  summarise(
    LN_IC50 = median(LN_IC50, na.rm = TRUE),
    .groups = "drop"
  )

ic50_mat <- ic50_summary %>%
  pivot_wider(
    names_from  = CELL_LINE_NAME,
    values_from = LN_IC50
  ) %>%
  column_to_rownames("DRUG_NAME") %>%
  as.matrix()

# 3. Rank the IC50 matrix
ic50_mat_ranked <- t(apply(ic50_mat, 1, rank, ties.method = "average"))

common_cells <- intersect(colnames(expr_mat_ranked), colnames(ic50_mat_ranked))
expr_final <- expr_mat_ranked[, common_cells]
ic50_final <- ic50_mat_ranked[, common_cells]

# 4. Correlation
spearman_cor_matrix <- cor(t(expr_final), t(ic50_final), method = "spearman", use = "pairwise.complete.obs")


library(reshape2)
cor_df <- melt(spearman_cor_matrix, varnames = c("gene", "drug"), value.name = "spearman_rho")

top_hits <- cor_df %>%
  filter(abs(spearman_rho) > 0.5) %>%
  arrange(desc(abs(spearman_rho)))

head(top_hits)
