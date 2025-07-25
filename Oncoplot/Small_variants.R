library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

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


#5. Define colors for mutation types
mut_colors <- c(
  "None" = "white",
  "Missense Mutation" = "#1f78b4",
  "Nonsense Mutation" = "#e31a6c",
  "Frame Shift" = "#33a02c",
  "In Frame" = "#b2df8a",
  "Splice Site" = "#f9df9b"
)

#5. Make sure mutation types in matrix are factor with these levels
mut_matrix[] <- lapply(mut_matrix, function(x) factor(x, levels = names(mut_colors)))

#6. Extract sample order from the existing ht_list
sample_order <- rownames(ht_list@ht_list[[1]]@matrix)

#7. Reorder mut_matrix to match that sample order
mut_matrix <- mut_matrix[sample_order, , drop = FALSE]

#8. Add a heatmap
ht_list <- ht_list + Heatmap(as.matrix(mut_matrix),
                             name = "Mutation Type",
                             col = mut_colors,
                             cluster_rows = FALSE,
                             cluster_columns = FALSE,
                             show_row_names = TRUE,
                             show_column_names = TRUE,
                             column_names_gp    = gpar(fontsize = 7),
                             na_col = "white",
                             rect_gp = gpar(col = "black", lwd = 0.5),
                             width              = unit(5, "cm"),
                             height             = unit(49 * 0.15, "cm"),
                             heatmap_legend_param = list(title = "Mutation Type", at = names(mut_colors), labels = names(mut_colors),
                                                         title_gp     = gpar(fontsize = 14),    
                                                         labels_gp    = gpar(fontsize = 14)),
                             
)


