library(readxl)
library(ggplot2)
library(patchwork)
library(dplyr)



data1 <- read_excel("./data/GDSC1_fitted_dose_response_27Oct23.xlsx") 

data2 <- read_excel("./data/GDSC2_fitted_dose_response_27Oct23.xlsx") 

# dimentions totale
dim(data1)
dim(data2)


# nombre de cell lines différentes
dim(table(data1$CELL_LINE_NAME))
dim(table(data2$CELL_LINE_NAME))


# distribution des categories TGCA (données groupées par cell line)
grouped_data1 <- data1 %>%
  group_by(CELL_LINE_NAME) %>%
  summarise(Count = n(),TCGA_DESC = paste(unique(TCGA_DESC), collapse = ", "))

grouped_data2 <- data2 %>%
  group_by(CELL_LINE_NAME) %>%
  summarise(Count = n(),TCGA_DESC = paste(unique(TCGA_DESC), collapse = ", "))

plot1 <- ggplot(grouped_data1, aes(x = TCGA_DESC, fill = TCGA_DESC)) +
  geom_bar() +
  theme_minimal() +
  labs(title = "Histogram of cell lines per TCGA_DESC Categories in GDSC1 (grouped by cell line)",
       x = "TCGA_DESC",
       y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none")

plot2 <- ggplot(grouped_data2, aes(x = TCGA_DESC, fill = TCGA_DESC)) +
  geom_bar() +
  theme_minimal() +
  labs(title = "Histogram of cell lines per TCGA_DESC Categories in GDSC2 (grouped by cell line)",
       x = "TCGA_DESC",
       y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none")

plot1 / plot2 + plot_layout(guides = "collect") & theme(legend.position = "none")


# comparaison des bdd par categorie TGCA (données groupées par cell line)

grouped_data1$Source <- "GDSC 1"
grouped_data2$Source <- "GDSC 2"
combined_data <- bind_rows(grouped_data1, grouped_data2)

all_data <- combined_data %>%
  group_by(Source,TCGA_DESC) %>%
  summarise(count = n(), .groups = 'drop')

total_sum <- sum(all_data$count)
all_data <- all_data %>%
  group_by(TCGA_DESC) %>%
  mutate(percentage = count / total_sum * 100)

ggplot(all_data, aes(x = TCGA_DESC, y = percentage, fill = Source)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Percentage distribution of GDSC 1 and GDSC 2 by TCGA Category (grouped by cell line)",
       x = "TCGA Category", 
       y = "Percentage",
       fill = "Database (GDSC)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Pour chaque drug, affiche le % de cell lines testé
combined_data <- bind_rows(data1, data2)
total_cell_lines <- length(unique(combined_data$CELL_LINE_NAME))

drug_cell_line_counts <- combined_data %>%
  group_by(DRUG_NAME) %>%
  summarise(test_count = n_distinct(CELL_LINE_NAME), .groups = 'drop')


drug_cell_line_counts <- drug_cell_line_counts %>%
  mutate(percentage = (test_count / total_cell_lines) * 100)


# Plot the results
ggplot(drug_cell_line_counts, aes(x = DRUG_NAME, y = percentage)) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "Percentage of Cell Lines Tested for Each Drug",
       x = "Drug name",
       y = "Percentage of Cell Lines Tested (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),legend.position = "none")


# comparaison des bdd par metabolic pathway (données groupées par cell line)
data1$Source <- "GDSC 1"
data2$Source <- "GDSC 2"
combined_pathway <- bind_rows(data1, data2)


all_data <- combined_pathway %>%
  group_by(Source,PATHWAY_NAME) %>%
  summarise(count = n(), .groups = 'drop') 


ggplot(all_data, aes(x = reorder(PATHWAY_NAME, -count), y = count, fill = Source)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Percentage distribution of GDSC 1 and GDSC 2 by metabolic pathway Category",
       x = "PATHWAY NAME", 
       y = "Count",
       fill = "Database (GDSC)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
