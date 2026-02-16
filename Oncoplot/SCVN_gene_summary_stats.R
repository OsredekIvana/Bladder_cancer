library(dplyr)

# 1. Count number of samples per gene with any alteration
gene_alteration_counts <- clean_scna %>%
  filter(!is.na(cnStatus) & cnStatus != "None") %>%
  distinct(Tumor_Sample_Barcode, gene) %>%
  count(gene, name = "n_samples_with_alteration")

# 2. Total number of unique samples
n_total_samples <- 23

# 3. Add percentage
gene_alteration_counts <- gene_alteration_counts %>%
  mutate(percent_altered = round(100 * n_samples_with_alteration / n_total_samples, 1)) %>%
  arrange(desc(percent_altered))

# 4. Count mutation types per gene
mutation_type_summary <- clean_scna %>%
  filter(!is.na(cnStatus) & cnStatus != "None") %>%
  count(gene, cnStatus) %>%
  group_by(gene) %>%
  mutate(total = sum(n),
         percent = round(100 * n / total, 1)) %>%
  ungroup()

# 5. Top genes summary with most frequent mutation type
top_gene_summary <- gene_alteration_counts %>%
  left_join(
    mutation_type_summary %>%
      group_by(gene) %>%
      slice_max(order_by = n, n = 1) %>%
      select(gene, most_common_type = cnStatus, n_type = n, percent_type = percent),
    by = "gene"
  )

# Print top 10
top_gene_summary %>% head(10)

# Load the library
library(writexl)

# Write the Excel file
write_xlsx(top_gene_summary, "/Users/osredek/Desktop/EMBL/PROJECT_Bladder_cancer/Oncoplots/Supplementary_Table_X.xlsx")



write.csv(top_gene_summary, "/Users/osredek/Desktop/EMBL/PROJECT_Bladder_cancer/Oncoplots/Supplementary_Table_X.csv", row.names = FALSE)

