library(ComplexHeatmap)
library(dplyr)
library(stringr)
library(tidyr)
library(circlize)
library(pheatmap)
library(data.table)
library(textshape)
library(circlize)


#1. Load data
df_clean <- read_excel("df_clean.xlsx")

#2. Convert SampleID column to rownames
df_clean <- as.data.frame(df_clean)
rownames(df_clean) <- df_clean$SampleID
df_clean$SampleID <- NULL

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

#3. Initialize the oncoplot 
ht_list = Heatmap(
  as.matrix(df_clean["Grade"]),
  name = "Grade",
  width = unit(0.4, "cm"),
  col = col$grade_colors,
  show_row_names = FALSE,
  border = TRUE,
  height = unit(18, "cm"),
  rect_gp = gpar(col = "white", lwd = 1),
  heatmap_legend_param = list(labels_gp = gpar(fontsize = 14),
                              title_gp = gpar(fontsize = 14))
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
                              title_gp = gpar(fontsize = 14)) 
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

#4. Add amplicon high/low status (manually currated with the cutoff <25)
amplicon_high <- c("B42","B123","B4","B5")
amplicon_low <- c("B134","B94","B125","B157","B187","B175","B12","B22",
                  "B32","B24","B154","B178","B67","B156","B39","B60","B87","B74")

amplicon_status <- setNames(rep("High", length(amplicon_high)), amplicon_high)
amplicon_status <- c(amplicon_status, setNames(rep("Low", length(amplicon_low)), amplicon_low))
status_vector <- amplicon_status[df_clean$rn]  
col$amplicon_status <- c("High" = "#E69F00", "Low" = "#56B4E9") 

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


#5. Read in the smoking data 
smoking_df <- read.csv("smoking_metadata.csv", stringsAsFactors = FALSE)

#6. Clean up the sample names 
smoking_df <- smoking_df %>%
  mutate(
    Tumor_Sample_Barcode = Tumor_Sample_Barcode %>%
      str_remove(regex("tumor", ignore_case = TRUE)) %>%
      str_replace_all("^B060$", "B60")
  )

#6. Merge into the clean_df
df_clean <- df_clean %>%
  tibble::rownames_to_column("Tumor_Sample_Barcode") %>%
  left_join(smoking_df, by = "Tumor_Sample_Barcode") %>%
  tibble::column_to_rownames("Tumor_Sample_Barcode")


#7.  Define colors for smoking status
smoke_cols <- c("YES" = "#D95F02", "NO" = "#7570B3")

#8. Append the track
ht_list <- ht_list +
  Heatmap(
    df_clean$Smoking_status.y,
    name               = "Smoking",
    col                = smoke_cols,
    width              = unit(0.4, "cm"),
    height             = unit(18, "cm"),      
    show_row_names     = TRUE,
    show_column_names  = TRUE,
    rect_gp            = gpar(col = "white", lwd = 1),
    heatmap_legend_param = list(
      labels_gp = gpar(fontsize = 14),
      title_gp  = gpar(fontsize = 14)
    )
  )

#9. Append the track for the fusions 
ht_list <- ht_list +
  Heatmap(
    df_clean$ZYG11B_RAF1,
    name               = "ZYG11B-RAF1",
    col                = smoke_cols,
    width              = unit(0.4, "cm"),
    height             = unit(18, "cm"),  
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
    col                = smoke_cols,
    width              = unit(0.4, "cm"),
    height             = unit(18, "cm"),      
    show_row_names     = TRUE,
    show_column_names  = TRUE,
    rect_gp            = gpar(col = "white", lwd = 1),
    heatmap_legend_param = list(
      labels_gp = gpar(fontsize = 14),
      title_gp  = gpar(fontsize = 14)
    )
  )


col_fun_dummy = circlize::colorRamp2(c(min(df_clean[["LINE1 insertions"]], na.rm = T), 49), c("white", "red"))

ht_list = ht_list + Heatmap(
  as.matrix(df_clean["LINE1 insertions"]),
  name = "LINE1 insertions",
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
    if (!is.na(df_clean[["LINE1 insertions"]][i])) {
      grid.text(df_clean[["LINE1 insertions"]][i], x, y, , rot = 90, gp = gpar(fontsize = 7, col = "black"))
    }
    # If it's NA, we don't draw text, and na_col makes the tile white
  },
  heatmap_legend_param = list(labels_gp = gpar(fontsize = 14),
                              title_gp = gpar(fontsize = 14))
)



#Add small mutations to the plot: See Small_variants.R
#Add SCNV mutations to the plot: See SCNV.R
