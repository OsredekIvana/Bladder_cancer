library(dplyr)
library(writexl)


#1. Read in the file 
raw_lines <- readLines("scna.gene.summary")

#2. Split by tabs or multiple spaces, remove empty strings
split_lines <- strsplit(raw_lines, "\\s+")

#3. Keep only rows with 6 fields
cleaned_lines <- split_lines[sapply(split_lines, length) == 6]

#4. Convert to data frame
scna <- as.data.frame(do.call(rbind, cleaned_lines), stringsAsFactors = FALSE)

#5. Assign column names
colnames(scna) <- c("sample", "cnStatus", "gene", "chr", "start", "end")

#6. Convert numeric columns
scna$start <- as.numeric(scna$start)
scna$end <- as.numeric(scna$end)

#7. Clean the dataframe
clean_scna <- scna %>%
  mutate(
    Tumor_Sample_Barcode = sample %>%
      str_remove(regex("tumor", ignore_case = TRUE)) %>%
      str_replace_all("^B060$", "B60"),
    
    cnStatus = recode(cnStatus,
                      "subclonalAmp"  = "Subclonal Amplification",
                      "cnLOH"         = "CN LOH",
                      "subclonalLOH"  = "Subclonal LOH",
                      "subclonalDel"  = "Subclonal Deletion"
    )
  ) %>%
  select(Tumor_Sample_Barcode, cnStatus, gene)

clean_scna <- clean_scna %>%
  filter(!gene %in% c("SPAN1", "NF1", "KMT2C", "ATM"))


#8. Count number of samples per gene with any alteration
gene_alteration_counts <- clean_scna %>%
  filter(!is.na(cnStatus) & cnStatus != "None") %>%
  distinct(Tumor_Sample_Barcode, gene) %>%
  count(gene, name = "n_samples_with_alteration")

#9. Total number of unique samples
n_total_samples <- 23

#10. Add percentage
gene_alteration_counts <- gene_alteration_counts %>%
  mutate(percent_altered = round(100 * n_samples_with_alteration / n_total_samples, 1)) %>%
  arrange(desc(percent_altered))

#11. Count mutation types per gene
mutation_type_summary <- clean_scna %>%
  filter(!is.na(cnStatus) & cnStatus != "None") %>%
  count(gene, cnStatus) %>%
  group_by(gene) %>%
  mutate(total = sum(n),
         percent = round(100 * n / total, 1)) %>%
  ungroup()

#12. Top genes summary with most frequent mutation type
top_gene_summary <- gene_alteration_counts %>%
  left_join(
    mutation_type_summary %>%
      group_by(gene) %>%
      slice_max(order_by = n, n = 1) %>%
      select(gene, most_common_type = cnStatus, n_type = n, percent_type = percent),
    by = "gene"
  )


#13. Write the Excel file
write.csv(top_gene_summary, "Supplementary_Table_D.csv", row.names = FALSE)

