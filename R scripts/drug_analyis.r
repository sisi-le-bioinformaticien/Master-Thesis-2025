library(readxl)
library(ggplot2)
library(patchwork)
library(tidyr)
library(dplyr)
library(DESeq2)
library(pheatmap)

MUT_DATA <- read.csv("./data/mutations_all.csv", sep=",") 
SMALL_DATA <- read.csv("./data/rnaseq_merged_20250117_filtered.csv",sep=",") 
data1 <- read_excel("./data/GDSC1_fitted_dose_response_27Oct23.xlsx") 
data2 <- read_excel("./data/GDSC2_fitted_dose_response_27Oct23.xlsx") 
combined_data <- bind_rows(data1, data2)

GENE_ID <- "ENSG00000173638" #SLC19A1/hRFC
DRUG_ENTRY <- "Methotrexate"

GENE_NAME <- "ABCC1"
DRUG_ENTRY <- "Methotrexate"
#GENE_ID <- "ENSG00000177606" 

GENE_ID <- MUT_DATA %>%
  filter(gene_symbol==GENE_NAME) %>%
  select(ensembl_gene_id) %>%
  pull(ensembl_gene_id) %>%
  .[1]


gdsc_cell_lines <- as.vector(unique(combined_data$CELL_LINE_NAME))


missing_cell_lines <- setdiff(gdsc_cell_lines, MUT_DATA$model_name)


MUT_DATAFRAME <- MUT_DATA %>%
  select(model_name, gene_symbol,ensembl_gene_id, protein_mutation, coding) %>%
  filter(model_name %in% gdsc_cell_lines, coding=="t") 


MUT_GENE <- MUT_DATAFRAME %>%
  select(model_name, gene_symbol,ensembl_gene_id, protein_mutation, coding) %>%
  filter(grepl(GENE_ID, ensembl_gene_id))

print(paste("number of cell lines mutatated in ", GENE_ID , ": " , length(table(unique(MUT_GENE$model_name)))))
print(paste("number of unique mutations in ", GENE_ID , ": " , length(table(unique(MUT_GENE$protein_mutation)))))


# PLOTTING each CELL LINE IC50 value for ENSG00000173638 gene and plot mutations vs wildtype

drug_ic50 <- combined_data %>%
  filter(DRUG_NAME == DRUG_ENTRY) %>%
  select(CELL_LINE_NAME, LN_IC50) 

drug_ic50 <- drug_ic50 %>%
  mutate(mutation_status = ifelse(CELL_LINE_NAME %in% MUT_GENE$model_name, "Mutated", "WT"))

ggplot(drug_ic50, aes(x = mutation_status, y = LN_IC50, color = mutation_status)) +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_boxplot(outlier.shape = NA, fill = NA) +
  theme_minimal() +
  labs(
    title = paste(DRUG_ENTRY, "sensitivity by", GENE_NAME, "Mutations"),
    x = "Mutation Status",
    y = "log(IC50)"
  )

# PLOTTING each cell line IC50 value in relation ENSG00000173638 with gene expression

data_cisplatin <- combined_data %>%
  filter(grepl(DRUG_ENTRY, DRUG_NAME)) %>%
  select(CELL_LINE_NAME,Z_SCORE)

data_cisplatin <- data_cisplatin %>%
  mutate(condition = case_when(
    Z_SCORE > 2  ~ "resistant",
    Z_SCORE < -2 ~ "sensitive",
    TRUE         ~ "other"
  ))


data_drug <- combined_data %>%
  filter(grepl(DRUG_ENTRY, DRUG_NAME)) %>%
  select(CELL_LINE_NAME, LN_IC50)

gdsc_cell_lines <- data_drug %>%
  pull(CELL_LINE_NAME)


DRUG_SPECIFIC_DATA <- SMALL_DATA %>%
  filter(model_name %in% gdsc_cell_lines) %>%
  filter(grepl(GENE_ID, ensembl_gene_id))

Gene_plot_data <- DRUG_SPECIFIC_DATA %>%
  select(model_name, expression = htseq_read_count) %>%
  inner_join(data_drug, by = c("model_name" = "CELL_LINE_NAME"))

Gene_plot_data <- Gene_plot_data %>%
  mutate(log_expr = log2(expression + 1)) %>%
  left_join(data_cisplatin, by = c("model_name" = "CELL_LINE_NAME"))


ggplot(Gene_plot_data, aes(x = log_expr, y = LN_IC50, color = condition)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  scale_color_manual(values = c(
    "resistant" = "firebrick",
    "sensitive" = "dodgerblue",
    "other" = "gray70"
  )) +
  theme_minimal() +
  labs(
    title = paste(DRUG_ENTRY, "Sensitivity vs", GENE_NAME, "Expression"),
    x = "log2(Expression + 1)",
    y = "log(IC50)",
    color = "Condition"
  )

filtred_data <- Gene_plot_data %>%
  filter(log_expr<8)


