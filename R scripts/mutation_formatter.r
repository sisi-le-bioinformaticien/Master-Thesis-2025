library(readxl)
library(ggplot2)
library(patchwork)
library(tidyr)
library(dplyr)
library(DESeq2)
library(pheatmap)

MUT_DATA <- read.csv("./data/mutations_all.csv", sep=",") 
data1 <- read_excel("./data/GDSC1_fitted_dose_response_27Oct23.xlsx") 
data2 <- read_excel("./data/GDSC2_fitted_dose_response_27Oct23.xlsx") 
combined_data <- bind_rows(data1, data2)

gdsc_cell_lines <- as.vector(unique(combined_data$CELL_LINE_NAME))


missing_cell_lines <- setdiff(gdsc_cell_lines, MUT_DATA$model_name)


MUT_DATA_GROUPED <- MUT_DATA %>%
  select(model_name, gene_symbol,ensembl_gene_id, protein_mutation, coding) %>%
  filter(!is.na(protein_mutation),coding=="t",model_name %in% gdsc_cell_lines) 

gene_mut_counts <- MUT_DATA_GROUPED %>%
  group_by(gene_symbol) %>%
  summarise(num_cell_lines = n_distinct(model_name)) %>%
  arrange(desc(num_cell_lines)) %>%
  filter(num_cell_lines > 5)

FINAL_MUT <- MUT_DATA_GROUPED %>%
  select(model_name, gene_symbol, ensembl_gene_id, protein_mutation) %>%
  filter(gene_symbol %in% gene_mut_counts$gene_symbol)

DATAFRAME <- FINAL_MUT %>%
  group_by(model_name, gene_symbol, ensembl_gene_id) %>%
  summarise(
    mutations = paste(protein_mutation, collapse = ";"),
    count_mutations = n(),
    .groups = 'drop'
  )

COUNT_GENES_PER_CELL_LINE <- DATAFRAME %>%
  group_by(model_name) %>%
  summarise(
    gene_count = n(),
    .groups = 'drop'
  ) %>%
  arrange(desc(gene_count))

barplot(COUNT_GENES_PER_CELL_LINE$gene_count, main = "Count of mutated genes per cell line", ylab ="gene count", xlab="cell lines")


HIGHMUT <- DATAFRAME %>%
  select(model_name, gene_symbol,ensembl_gene_id, mutations, count_mutations)%>%
  filter(count_mutations>1) %>%
  arrange(desc(count_mutations))



gene_size <- read.csv("./data/gene_size_basic.csv")

gene_size <- gene_size %>%
  mutate(gene_id = sub("\\..*", "", gene_id))

NORMALIZED_HIGH_MUT <- HIGHMUT %>%
  left_join(gene_size, by = c("ensembl_gene_id" = "gene_id")) %>%
  mutate(
    normalized_mutation_rate = (count_mutations / size) * 1000
  ) %>%
  arrange(desc(normalized_mutation_rate))

