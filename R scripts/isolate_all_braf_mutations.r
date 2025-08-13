library(dplyr)


MUT_DATA <- read.csv("./data/mutations_all.csv", sep=",") 

BRAF_profiles <- MUT_DATA %>%
  filter(gene_symbol == "BRAF") %>%
  mutate(protein_mutation = ifelse(
    protein_mutation %in% c("-", "p.?", "", NA), 
    NA, 
    protein_mutation
  )) %>%
  group_by(model_name) %>%
  summarise(BRAF_mutations = paste(sort(na.omit(unique(protein_mutation))), collapse = "; ")) %>%
  mutate(BRAF_mutations = ifelse(BRAF_mutations == "", "Wildtype", BRAF_mutations)) %>%
  ungroup()

BRAF_profiles_summary <- BRAF_profiles %>%
  group_by(BRAF_mutations) %>%
  summarise(count = n()) %>%
  arrange(desc(count))                

write.csv(BRAF_profiles, 
          file = "BRAF_cell_line_mut_index.csv", 
          row.names = FALSE)
write.csv(BRAF_profiles_summary, 
          file = "BRAF_mutations_count.csv", 
          row.names = FALSE)

cat("Total unique BRAF mutations:", nrow(BRAF_mutations), "\n")
cat("Top 10 most frequent BRAF mutations:\n\n")

top_10 <- head(BRAF_mutations, 10)
for(i in 1:nrow(top_10)) {
  cat(paste0(i, ". ", top_10$protein_mutation[i], 
             " -> ", top_10$count[i], " occurrences\n"))
}
