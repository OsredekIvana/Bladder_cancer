library(dplyr)
library(tidyr)
library(writexl)


#1. Read data 
df_cleaned <- read.csv("cleaned_maf_df_for_ComplexHeatmap.csv", stringsAsFactors = FALSE)

df_cleaned <- df_cleaned %>%
  mutate(
    Tumor_Sample_Barcode = Tumor_Sample_Barcode %>%
      str_remove(regex("tumor", ignore_case = TRUE)) %>%
      str_replace_all("^B060$", "B60")
  )


#2. Select relevant columns and keep one mutation type per sample/gene
gene_sample_mut <- df_cleaned %>%
  select(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification) %>%
  mutate(Variant_Classification = recode(Variant_Classification,
                                         "Missense_Mutation" = "Missense Mutation",
                                         "Nonsense_Mutation" = "Nonsense Mutation",
                                         "Frame_Shift"       = "Frame Shift",
                                         "In_Frame"          = "In Frame",
                                         "Splice_Site"          = "Splice Site"
  )) %>%
  distinct() %>%
  group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%
  ungroup() 

#3. Pivot to wide
mut_matrix <- gene_sample_mut %>%
  pivot_wider(
    names_from = Hugo_Symbol,
    values_from = Variant_Classification,
    values_fill = NA
  ) %>%
  column_to_rownames("Tumor_Sample_Barcode")




#4. Count number of samples per gene with any alteration
gene_alteration_counts_SNV <- gene_sample_mut %>%
  filter(!is.na(Variant_Classification) & Variant_Classification != "None") %>%
  distinct(Tumor_Sample_Barcode, Hugo_Symbol) %>%
  count(Hugo_Symbol, name = "n_samples_with_alteration")

#5. Total number of unique samples
n_total_samples <- 23

#6. Add percentage
gene_alteration_counts_SNV <- gene_alteration_counts_SNV %>%
  mutate(percent_altered = round(100 * n_samples_with_alteration / n_total_samples, 1)) %>%
  arrange(desc(percent_altered))


#7. Count mutation types per gene
mutation_type_summary <- gene_sample_mut %>%
  filter(!is.na(Variant_Classification) & Variant_Classification != "None") %>%
  count(Hugo_Symbol, Variant_Classification) %>%
  group_by(Hugo_Symbol) %>%
  mutate(total = sum(n),
         percent = round(100 * n / total, 1)) %>%
  ungroup()

#8. Top genes summary with most frequent mutation type
top_gene_summary <- gene_alteration_counts_SNV %>%
  left_join(
    mutation_type_summary %>%
      group_by(Hugo_Symbol) %>%
      slice_max(order_by = n, n = 1) %>%
      select(Hugo_Symbol, most_common_type = Variant_Classification, n_type = n, percent_type = percent),
    by = "Hugo_Symbol"
  )


#9. Write the Excel file
write.csv(top_gene_summary, "Supplementary_Table_C.csv", row.names = FALSE)

