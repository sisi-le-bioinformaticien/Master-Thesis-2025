library(tidyverse)


braf_mutations <- read.csv("BRAF_cell_line_mut_index2.csv") %>%
  mutate(model_name = toupper(model_name))

data1 <- readxl::read_excel("./data/GDSC1_fitted_dose_response_27Oct23.xlsx")
data2 <- readxl::read_excel("./data/GDSC2_fitted_dose_response_27Oct23.xlsx")
combined_data <- bind_rows(data1, data2) %>%
  mutate(CELL_LINE_NAME = toupper(CELL_LINE_NAME)) 

DRUG_ENTRY <- "Dabrafenib"  

# Create mutation categories
mut_categories <- braf_mutations %>%
  mutate(
    mutation_group = case_when(
      BRAF_mutations == "Wildtype" ~ "Wildtype",
      grepl("V600E", BRAF_mutations) ~ "V600E",
      TRUE ~ "Others"
    )
  ) %>%
  select(model_name, mutation_group)

# Get IC50 values for the specified drug
drug_ic50 <- combined_data %>%
  filter(DRUG_NAME == DRUG_ENTRY) %>%
  select(CELL_LINE_NAME, LN_IC50) %>%
  distinct(CELL_LINE_NAME, .keep_all = TRUE)  # Remove duplicates

# Merge mutation data with IC50 values
plot_data <- drug_ic50 %>%
  left_join(mut_categories, by = c("CELL_LINE_NAME" = "model_name")) %>%
  filter(!is.na(mutation_group))  # Remove cell lines without mutation data

# Create boxplot
ggplot(plot_data, aes(x = mutation_group, y = LN_IC50, fill = mutation_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  scale_fill_manual(values = c("Wildtype" = "grey80", "V600E" = "#FF9999", "Others" = "#66B2FF")) +
  theme_minimal() +
  labs(
    title = paste(DRUG_ENTRY, "Sensitivity by BRAF Mutation Status"),
    subtitle = paste("Cell lines:", nrow(plot_data)),
    x = "BRAF Mutation Category",
    y = "log(IC50)",
    fill = "Mutation Group"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

write.csv(plot_data, "BRAF_mutations_with_IC50.csv", row.names = FALSE)