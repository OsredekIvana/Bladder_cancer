library(dplyr)
library(tidyr)


#Read in the file 
raw_lines <- readLines("scna.gene.summary")

# Split by tabs or multiple spaces, remove empty strings
split_lines <- strsplit(raw_lines, "\\s+")

# Keep only rows with 6 fields
cleaned_lines <- split_lines[sapply(split_lines, length) == 6]

# Convert to data frame
scna <- as.data.frame(do.call(rbind, cleaned_lines), stringsAsFactors = FALSE)

# Assign column names manually
colnames(scna) <- c("sample", "cnStatus", "gene", "chr", "start", "end")

# Convert numeric columns
scna$start <- as.numeric(scna$start)
scna$end <- as.numeric(scna$end)

####
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

filtered_scna <- clean_scna %>%
  filter(!gene %in% c("SPAN1", "NF1", "KMT2C", "ATM"))


#Include in the heatmap
# Pivot filtered SCNA to wide format: rows = samples, columns = genes
scna_matrix <- filtered_scna %>%
  select(Tumor_Sample_Barcode, gene, cnStatus) %>%
  distinct() %>%
  pivot_wider(
    names_from = gene,
    values_from = cnStatus,
    values_fill = NA
  ) %>%
  column_to_rownames("Tumor_Sample_Barcode")

scna_colors <- c(
  "Amplification" = "#e41a1c",
  "Aubclonal Amplification" = "#ff7f00",
  "CN LOH"        = "#377eb8",
  "Subclonal LOH" = "#984ea3",
  "Deletion"     = "#4daf4a",
  "Subclonal Del" = "#a65628",
  "None"         = "white"  # For empty/NA cells
)

# Convert values to factors
scna_matrix[] <- lapply(scna_matrix, function(x) factor(x, levels = names(scna_colors)))

# Ensure SCNA matrix matches the sample order from ht_list
sample_order <- rownames(ht_list@ht_list[[1]]@matrix)
scna_matrix <- scna_matrix[sample_order, , drop = FALSE]

ht_list <- ht_list + Heatmap(as.matrix(scna_matrix),
                             name = "SCNA Status",
                             col = scna_colors,
                             cluster_rows = FALSE,
                             cluster_columns = FALSE,
                             show_row_names = FALSE,
                             show_column_names = TRUE,
                             column_names_gp = gpar(fontsize = 7),
                             na_col = "white",
                             rect_gp = gpar(col = "black", lwd = 0.5),
                             width = unit(10, "cm"),
                             height = unit(49 * 0.15, "cm"),
                             heatmap_legend_param = list(
                               title = "SCNA Status",
                               at = names(scna_colors),
                               labels = names(scna_colors),
                               title_gp = gpar(fontsize = 14),
                               labels_gp = gpar(fontsize = 14)
                             )
)
`