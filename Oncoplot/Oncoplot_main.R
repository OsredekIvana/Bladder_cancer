library(readxl)
library(ComplexHeatmap)
library(dplyr)
library(stringr)
library(tidyr)
library(circlize)
library(pheatmap)
library(data.table)
library(maftools)
library(grid)

# Read in the Excel file
df_clean_in <- read_excel("df_clean.xlsx")

# Convert SampleID column to rownames
df_clean <- as.data.frame(df_clean_in)
rownames(df_clean) <- df_clean$SampleID
df_clean$SampleID <- NULL
df_clean_in <- df_clean

col <- list(
  sex_color = c("male" = "#56B4E9", "female" = "pink"),
  stage_colors = c(
    "Benign" = "coral1",
    "pTa" = "#E69F00",
    "pTis" = "#F0E442",
    "pT1" = "#CC79A7",
    "pT2" = "#999999"
  ),
  grade_colors = c(
    "Benign" = "coral1",
    "Low Grade" = "#41afaa",
    "High Grade" = "darkorchid"
  ),
  sequencing_colors = c("YES" = "dimgrey", "NO" = "azure2"),
  ecDNA_colors = c("NO" = "#949494", "YES" = "#78D3D3"),
  tp53_colors=c("None" = "white", "LOH" = "#D55E60", "Duplication" = "#CC79A7", "Deletion" = "black"),
  ZYG11B_RAF1=c("NO" = "#949994", "YES" = "#11D1f1"),
  FGFR3_TACC3=c("NO" = "#949911", "YES" = "#113155")
)

# Create the data frame
tp53_status <- data.frame(
  Tumor_Sample_Barcode = c("B4", "B5", "B12", "B22", "B24", "B32", "B39", "B42", "B60",
                           "B67", "B74", "B85", "B87", "B94", "B123", "B125", "B134", 
                           "B154", "B156", "B157", "B175", "B178", "B187"),
  TP53 = c("None","LOH","LOH","None","LOH","LOH", "LOH", "LOH","None","LOH","LOH","LOH","None","None","LOH", 
           "None","Duplication","LOH","None","None","None","None","Deletion"))


df_clean$rn <- rownames(df_clean)

df_clean <- inner_join(df_clean, tp53_status,
                       by = c("rn" = "Tumor_Sample_Barcode"))

rownames(df_clean) <- df_clean$rn

# L1 group definition:
# L1-low < 5, L1-intermediate 5–15, L1-high > 15
df_clean$L1_group <- dplyr::case_when(
  df_clean[["L1 insertions"]] < 5  ~ "L1-low",
  df_clean[["L1 insertions"]] > 15 ~ "L1-high",
  TRUE                             ~ "L1-intermediate"
)

df_clean$L1_group <- factor(
  df_clean$L1_group,
  levels = c("L1-low", "L1-intermediate", "L1-high")
)

# Reorder to make association visually obvious
df_clean <- df_clean[order(df_clean$Grade, df_clean$L1_group), ]

ht_list = Heatmap(
  as.matrix(df_clean["Grade"]),
  name = "Grade",
  # row_order = row_order_grade,
  width = unit(0.4, "cm"),
  col = col$grade_colors,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_rot = 180,
  row_names_centered = T,
  border = TRUE,
  height = unit(18, "cm"),
  rect_gp = gpar(col = "white", lwd = 1),
  heatmap_legend_param = list(labels_gp = gpar(fontsize = 14),
                              title_gp = gpar(fontsize = 14))
)

col_fun_dummy = circlize::colorRamp2(c(min(df_clean[["L1 insertions"]], na.rm = T), 49), c("white", "red"))


ht_list = ht_list + Heatmap(
  as.matrix(df_clean["L1 insertions"]),
  name = "L1 insertions",
  show_row_names = TRUE,
  col = col_fun_dummy,
  width = unit(0.4, "cm"),
  row_names_rot = 180,
  row_names_centered = T,
  border = TRUE,
  height = unit(18, "cm"),
  rect_gp = gpar(col = "white", lwd = 1),
  cell_fun = function(j, i, x, y, width, height, fill) {
    # Only draw text if the value is NOT NA
    if (!is.na(df_clean[["L1 insertions"]][i])) {
      grid.text(df_clean[["L1 insertions"]][i], x, y, , rot = 90, gp = gpar(fontsize = 7, col = "black"))
    }
    # If it's NA, don't draw text, and na_col makes the tile white
  },
  heatmap_legend_param = list(labels_gp = gpar(fontsize = 14),
                              title_gp = gpar(fontsize = 14))
)



col$L1_group <- c("L1-low" = "#56B4E9", "L1-high" = "#E69F00", "L1-intermediate" = "#999999")

# Add L1 group bar
ht_list = ht_list + Heatmap(
  df_clean$L1_group,
  name = "L1 group",
  col = col$L1_group,
  width = unit(0.4, "cm"),
  border = TRUE,
  height = unit(18, "cm"),
  show_row_names = T,
  show_column_names = T,
  rect_gp = gpar(col = "white", lwd = 1),
  heatmap_legend_param = list(labels_gp = gpar(fontsize = 14),
                              title_gp  = gpar(fontsize = 14))
)


ht_list = ht_list + Heatmap(
  df_clean$stage,
  name = "Stage",
  width = unit(0.4, "cm"),
  col = col$stage_colors,
  rect_gp = gpar(col = "white", lwd = 1),
  border = TRUE,
  height = unit(18, "cm"),
  heatmap_legend_param = list(labels_gp = gpar(fontsize = 14),
                              title_gp = gpar(fontsize = 14))
)

ht_list = ht_list + Heatmap(
  df_clean$sex,
  name = "Sex",
  width = unit(0.4, "cm"),
  col = col$sex_color,
  border = TRUE,
  height = unit(18, "cm"),
  rect_gp = gpar(col = "white", lwd = 1),
  heatmap_legend_param = list(labels_gp = gpar(fontsize = 14),
                              title_gp = gpar(fontsize = 14)) # Adjust this value to increase or decrease text size)
)

ht_list = ht_list + Heatmap(
  df_clean$Cell.free.DnA,
  name = "cfDNA",
  width = unit(0.4, "cm"),
  col = col$sequencing_colors,
  border = TRUE,
  height = unit(18, "cm"),
  rect_gp = gpar(col = "white", lwd = 1),
  heatmap_legend_param = list(labels_gp = gpar(fontsize = 14),
                              title_gp = gpar(fontsize = 14))
)

ht_list = ht_list + Heatmap(
  df_clean$ONT.experiment,
  name = "ONT",
  width = unit(0.4, "cm"),
  col = col$sequencing_colors,
  border = TRUE,
  height = unit(18, "cm"),
  rect_gp = gpar(col = "white", lwd = 1),
  heatmap_legend_param = list(labels_gp = gpar(fontsize = 14),
                              title_gp = gpar(fontsize = 14))
)


ht_list = ht_list + Heatmap(
  df_clean$Bulk.RNA,
  name = "RNA-seq",
  width = unit(0.4, "cm"),
  col = col$sequencing_colors,
  border = TRUE,
  height = unit(18, "cm"),
  rect_gp = gpar(col = "white", lwd = 1),
  heatmap_legend_param = list(labels_gp = gpar(fontsize = 14),
                              title_gp = gpar(fontsize = 14))
)


ht_list = ht_list + Heatmap(
  as.matrix(df_clean["Visium"]),
  name = "Visium",
  width = unit(0.4, "cm"),
  height = unit(18, "cm"),
  col = col$sequencing_colors,
  show_row_names = FALSE,
  border = TRUE,
  row_names_rot = 180,
  row_names_centered = T,
  rect_gp = gpar(col = "white", lwd = 1),
  heatmap_legend_param = list(labels_gp = gpar(fontsize = 14),
                              title_gp = gpar(fontsize = 14))
)

ht_list = ht_list + Heatmap(
  df_clean$ecDNA,
  name = "ecDNA",
  width = unit(0.4, "cm"),
  col = col$ecDNA_colors,
  border = TRUE,
  height = unit(18, "cm"),
  rect_gp = gpar(col = "white", lwd = 1),
  heatmap_legend_param = list(labels_gp = gpar(fontsize = 14),
                              title_gp = gpar(fontsize = 14))
)


#Add amplicon high/low status 

amplicon_high <- c("B42","B123","B4","B5")
amplicon_low <- c("B134","B94","B125","B157","B187","B175","B12","B22",
                  "B32","B24","B154","B178","B67","B156","B39","B60","B87","B74","B85", "B154","B156")

amplicon_status <- setNames(rep("High", length(amplicon_high)), amplicon_high)
amplicon_status <- c(amplicon_status, setNames(rep("Low", length(amplicon_low)), amplicon_low))
status_vector <- amplicon_status[df_clean$rn]  # Adjust column name as needed
col$amplicon_status <- c("High" = "#E69F00", "Low" = "#56B4E9")  # red and blue

ht_list <- ht_list + 
  Heatmap(
    status_vector,
    name = "Amplicon",
    col = col$amplicon_status,
    width = unit(0.4, "cm"),
    border = TRUE,
    height = unit(18, "cm"),
    rect_gp = gpar(col = "white", lwd = 1),
    heatmap_legend_param = list(
      labels_gp = gpar(fontsize = 14),
      title_gp = gpar(fontsize = 14)
    )
  )


# Define colors for smoking status
smoke_cols <- c("YES" = "#D95F02", "NO" = "#7570B3")

# Append the track
ht_list <- ht_list +
  Heatmap(
    df_clean$Smoking_status,
    name               = "Smoking",
    col                = smoke_cols,
    width              = unit(0.4, "cm"),
    height             = unit(18, "cm"),       # match your other tracks
    show_row_names     = TRUE,
    show_column_names  = TRUE,
    rect_gp            = gpar(col = "white", lwd = 1),
    heatmap_legend_param = list(
      labels_gp = gpar(fontsize = 14),
      title_gp  = gpar(fontsize = 14)
    )
  )

# Append the track for the fusions 
ht_list <- ht_list +
  Heatmap(
    df_clean$ZYG11B_RAF1,
    name               = "ZYG11B-RAF1",
    col                = col$ZYG11B_RAF1,
    width              = unit(0.4, "cm"),
    height             = unit(18, "cm"),       # match your other tracks
    show_row_names     = TRUE,
    show_column_names  = TRUE,
    rect_gp            = gpar(col = "white", lwd = 1),
    heatmap_legend_param = list(
      labels_gp = gpar(fontsize = 14),
      title_gp  = gpar(fontsize = 14)
    )
  )

ht_list <- ht_list +
  Heatmap(
    df_clean$FGFR3_TACC3,
    name               = "FGFR3-TACC3",
    col                = col$FGFR3_TACC3,
    width              = unit(0.4, "cm"),
    height             = unit(18, "cm"),       # match your other tracks
    show_row_names     = TRUE,
    show_column_names  = TRUE,
    rect_gp            = gpar(col = "white", lwd = 1),
    heatmap_legend_param = list(
      labels_gp = gpar(fontsize = 14),
      title_gp  = gpar(fontsize = 14)
    )
  )


#Add the variants 

#Maf file preparation 

#Read maf 
maf_file = read.maf(maf="summary.maf",rmFlags=TRUE)

#Subset the CNV table to only those genes that are present in the oncoplot
gene_summary <- maf_file@gene.summary$Hugo_Symbol
cn.data_filtered_for_maf <- cn.data[cn.data$Gene %in% gene_summary, ]

str(gene_summary)

top_genes <- maf_file@gene.summary %>%
  arrange(desc(MutatedSamples)) %>%   # or desc(Total)
  slice_head(n = 20) %>%
  pull(Hugo_Symbol)

maf_df <- maf_file@data

top_mutations_df <- maf_df %>%
  filter(Hugo_Symbol %in% top_genes) %>%
  select(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification) 

#Merge following Variant_Classification rows into: 
#Frame_Shift_Del and Frame_Shift_Ins into Frame_Shift 
#In_Frame_Del and In_Frame_Ins into In_Frame
#If after merging, same gene for same sample results in duplicates, remove them

maf_df_cleaned <- top_mutations_df %>%
  mutate(
    Variant_Classification = case_when(
      Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins") ~ "Frame_Shift",
      Variant_Classification %in% c("In_Frame_Del", "In_Frame_Ins") ~ "In_Frame",
      TRUE ~ Variant_Classification
    )
  )

#Remove duplicates (same gene in same sample, post-merge)
maf_df_cleaned <- maf_df_cleaned %>%
  distinct(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, .keep_all = TRUE)


####Small variants#####

library(dplyr)
library(tidyr)
library(circlize)
library(textshape)

maf_df_cleaned <- read.csv(
  "cleaned_maf_df_for_ComplexHeatmap.csv",
  stringsAsFactors = FALSE
)

# 1. Select relevant columns and keep one mutation type per sample/gene (if multiple, pick first)
gene_sample_mut <- maf_df_cleaned %>%
  select(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification) %>%
  mutate(
    # clean sample names
    Tumor_Sample_Barcode = as.character(Tumor_Sample_Barcode),
    Tumor_Sample_Barcode = trimws(Tumor_Sample_Barcode),
    Tumor_Sample_Barcode = sub("tumor$", "", Tumor_Sample_Barcode),          # remove 'tumor'
    Tumor_Sample_Barcode = sub("^B0+([0-9]+)$", "B\\1", Tumor_Sample_Barcode), # remove extra zeros (B060 -> B60)
    
    # recode mutation types
    Variant_Classification = recode(Variant_Classification,
                                    "Missense_Mutation" = "Missense Mutation",
                                    "Nonsense_Mutation" = "Nonsense Mutation",
                                    "Frame_Shift"       = "Frame Shift",
                                    "In_Frame"          = "In Frame",
                                    "Splice_Site"       = "Splice Site"
    )
  ) %>%
  distinct() %>%
  group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%
  slice(1) %>%
  ungroup()

# 2. Pivot to wide: rows = samples, columns = genes, values = mutation types
mut_matrix <- gene_sample_mut %>%
  pivot_wider(
    names_from = Hugo_Symbol,
    values_from = Variant_Classification
  ) %>%
  column_to_rownames("Tumor_Sample_Barcode")



# Define colors for mutation types (add more as needed)
mut_colors <- c(
  "None" = "white",
  "Missense Mutation" = "#1f78b4",
  "Nonsense Mutation" = "#e31a6c",
  "Frame Shift" = "#33a02c",
  "In Frame" = "#b2df8a",
  "Splice Site" = "#f9df9b"
)

# Make sure mutation types in matrix are factor with these levels
mut_matrix[] <- lapply(mut_matrix, function(x) factor(x, levels = names(mut_colors)))

# Extract sample order from the existing ht_list
sample_order <- rownames(ht_list@ht_list[[1]]@matrix)

# Reorder mut_matrix to match that sample order (for some reason those samples were NAs, but I identified which is which)
mut_matrix <- mut_matrix[sample_order, , drop = FALSE]

rn <- rownames(mut_matrix)

rn[rn == "NA"]   <- "B178"
rn[rn == "NA.1"] <- "B74"

rownames(mut_matrix) <- rn


#Add a heatmap
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





# 1. Count number of samples per gene with any alteration
gene_alteration_counts_SNV <- gene_sample_mut %>%
  filter(!is.na(Variant_Classification) & Variant_Classification != "None") %>%
  distinct(Tumor_Sample_Barcode, Hugo_Symbol) %>%
  count(Hugo_Symbol, name = "n_samples_with_alteration")

# 2. Total number of unique samples
n_total_samples <- 23

# 3. Add percentage
gene_alteration_counts_SNV <- gene_alteration_counts_SNV %>%
  mutate(percent_altered = round(100 * n_samples_with_alteration / n_total_samples, 1)) %>%
  arrange(desc(percent_altered))


# 4. Count mutation types per gene
mutation_type_summary <- gene_sample_mut %>%
  filter(!is.na(Variant_Classification) & Variant_Classification != "None") %>%
  count(Hugo_Symbol, Variant_Classification) %>%
  group_by(Hugo_Symbol) %>%
  mutate(total = sum(n),
         percent = round(100 * n / total, 1)) %>%
  ungroup()

# 5. Top genes summary with most frequent mutation type
top_gene_summary <- gene_alteration_counts_SNV %>%
  left_join(
    mutation_type_summary %>%
      group_by(Hugo_Symbol) %>%
      slice_max(order_by = n, n = 1) %>%
      select(Hugo_Symbol, most_common_type = Variant_Classification, n_type = n, percent_type = percent),
    by = "Hugo_Symbol"
  )

# Print top 10
top_gene_summary %>% head(10)


# Load the library
library(writexl)

# Write the Excel file
write_xlsx(top_gene_summary, "Supplementary_Table_X1.xlsx")

write.csv(top_gene_summary, "Supplementary_Table_X1.csv", row.names = FALSE)


#######









