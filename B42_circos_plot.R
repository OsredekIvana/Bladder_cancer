library(dplyr)
library(circlize)
library(tidyr)
library(data.table)
library(GenomicRanges)
library(readxl)
library(fuzzyjoin)
library(readr)
library(ComplexHeatmap)

# 1. Read promoter data (manually generated)
B42_detailed_L1 <- read_excel("B42_detailed_L1.xlsx")

# 2. Prepare BED regions for promoters
bed_inGRCh38 <- B42_detailed_L1 %>%
  filter(L1 == "inGRCh38") %>%
  select(chr, promotor_start, promotor_end)

bed_notGRCh38 <- B42_detailed_L1 %>%
  filter(L1 %in% c("notInGRCh38", "somatic")) %>%
  filter(!is.na(promotor_start)) %>%
  select(chr, promotor_start, promotor_end)

# 3. Read Pileup for notInGRCh38 L1s
path_to_files <- "Pileup_files_notInGRCh38"
files <- list.files(path_to_files, pattern = "_assembly.pileup.bed.gz", full.names = TRUE)
read_pileup_file <- function(file) {
  if (file.info(file)$size == 0) {
    message(paste("Skipping empty file:", basename(file)))
    return(NULL)
  }
  tryCatch({
    tmp <- read_tsv(file, col_names = FALSE, comment = "#")
    if (ncol(tmp) == 18) {
      colnames(tmp) <- c("chrom", "start", "end", "modified_base", "score", "strand",
                         "start_incl", "end_incl", "color", "Nvalid_cov", "fraction_modified",
                         "Nmod", "Ncanonical", "Nother_mod", "Ndelete", "Nfail", "Ndiff", "Nnocall")
      tmp$file_name <- basename(file)
      return(tmp)
    } else {
      message(paste("Skipping file with unexpected column count:", basename(file)))
      return(NULL)
    }
  }, error = function(e) {
    message(paste("Failed to read file:", basename(file)))
    return(NULL)
  })
}

all_data <- lapply(files, read_pileup_file)
notInGRCh38 <- rbindlist(Filter(Negate(is.null), all_data), fill = TRUE)
notInGRCh38$file_name <- gsub("_assembly.pileup.bed.gz", "", notInGRCh38$file_name)

# 4. Read methylation data for inGRCh38 promoters and bodies
promoter_meth_inGRCh38_modkit <- fread("promoter_methylation_inGRCh38.bed", header = FALSE)
body_meth_inGRCh38_modkit <- fread("L1_promoter_methylation/body_methylation_inGRCh38.bed", header = FALSE)

# 5. Match promoter methylation positions to promoter regions using fuzzy join
matched_promoter <- fuzzy_inner_join(
  promoter_meth_inGRCh38_modkit,
  bed_inGRCh38,
  by = c("V1" = "chr", "V2" = "promotor_start", "V3" = "promotor_end"),
  match_fun = list(`==`, `>=`, `<=`)
)


promoter_meth_InGRCh38 <- matched_promoter %>%
  mutate(file_name = paste0(chr, "_", promotor_start, "_", promotor_end))

# 6. Prepare body regions for inGRCh38 L1s, handling chrX separately
bed_inGRCh38_body <- B42_detailed_L1 %>%
  filter(L1 == "inGRCh38", chr != "chrX") %>%
  select(chr, body_start, body_end)

bed_inGRCh38_body_chrX <- B42_detailed_L1 %>%
  filter(L1 == "inGRCh38", chr == "chrX") %>%
  select(chr, body_start, body_end)

matched_body <- fuzzy_inner_join(
  body_meth_inGRCh38_modkit,
  bed_inGRCh38_body,
  by = c("V1" = "chr"),
  match_fun = list(`==`)
)

# 7. Define promoter and body regions for notInGRCh38 L1s manually
promoter_not_inGRCh38 <- data.frame(
  file_name = c("chr11_59485054_59495054", "chr2_143248214_143258215", "chr2_190608314_190618314", 
                "chr3_123866865_123876865", "chr3_186652727_186662727", "chr5_108245797_108255797", 
                "chr5_110138384_110148384", "chr6_13185784_13195787", "chr6_19787888_19797891", 
                "chr9_112009267_112019273"),
  promoter_start = c(1123, 13098, 11018, 17559, 20267, 19391, 14406, 27767, 13924, 20842),
  promoter_end = c(1623, 13598, 11518, 18059, 20767, 19891, 14906, 28267, 14424, 21342)
)

body_not_inGRCh38 <- data.frame(
  file_name = c("chr11_59485054_59495054", "chr2_143248214_143258215", "chr2_190608314_190618314", 
                "chr3_123866865_123876865", "chr3_186652727_186662727", "chr5_108245797_108255797", 
                "chr5_110138384_110148384", "chr6_13185784_13195787", "chr6_19787888_19797891", 
                "chr9_112009267_112019273"),
  body_start = c(1624, 13599, 5506, 18060, 14752, 13876, 8891, 22253, 14425, 15328),
  body_end = c(7138, 19112, 11017, 23573, 20266, 19390, 14405, 27766, 19938, 20841)
)

# 8. Annotate methylation data with promoter/body coordinates for notInGRCh38 L1s
annotated_promoter <- merge(notInGRCh38, promoter_not_inGRCh38, by = "file_name", all.x = TRUE)
annotated_body <- merge(notInGRCh38, body_not_inGRCh38, by = "file_name", all.x = TRUE)

promoter_meth_notInGRCh38 <- annotated_promoter[start_incl >= promoter_start & start_incl <= promoter_end]
body_meth_notInGRCh38 <- annotated_body[start_incl >= body_start & start_incl <= body_end]

# 9. Calculate average methylation for bodies
avrg_body_meth_notInGRCh38 <- body_meth_notInGRCh38 %>%
  group_by(file_name) %>%
  summarise(avrg_body_meth_notInGRCh38 = mean(fraction_modified))

avrg_body_meth_inGRCh38 <- matched_body %>%
  group_by(V1) %>%
  summarise(avrg_body_meth_notInGRCh38 = mean(V11))

# 10. Handle chrX body methylation separately due to matching complexity
body_meth_inGRCh38_modkit_chrX117 <- body_meth_inGRCh38_modkit %>%
  filter(V1 == "chrX", V2 < 11900000) %>%
  summarise(avrg_body_meth_notInGRCh38 = mean(V11))

body_meth_inGRCh38_modkit_chrX119 <- body_meth_inGRCh38_modkit %>%
  filter(V1 == "chrX", V2 > 11900000) %>%
  summarise(avrg_body_meth_notInGRCh38 = mean(V11))


# 11. Read source/target partner L1 table
B42_file_source_target <- read.table(
  "L1_summary_B42.txt",
  header = FALSE,
  sep = "\t",
  col.names = c("chr", "pos", "type", "mei", "source1", "source2")
)

B42_partnered_L1s <- B42_file_source_target %>%
  filter(type == "partnered", mei == "L1")

# 12. Tables (manually currated) for methylation summary per source2 region (promoter and body)
source2_meth_promoter <- tibble::tibble(
  source2 = c(
    "chr11:24334007-24336603","chr11:59490054-59491029", "chr2:143253214-143253651",
    "chr2:143253214-143253674", "chr2:143253214-143253679", "chr2:143253215-143253655",
    "chr2:190613314-190613981", "chr22:28669325-28669371", "chr22:28669325-28669921",
    "chr22:28669325-28670131", "chr22:28669325-28670432", "chr22:28669541-28670131",
    "chr22:28669551-28670131", "chr3:123871865-123872261", "chr5:108250797-108251136",
    "chr5:110143384-110144564", "chr6:13190784-13190810", "chr6:13190784-13190820",
    "chr6:13190784-13191209", "chr6:13190784-13191211", "chr6:13190784-13191212",
    "chr6:13190787-13190817", "chr6:13190787-13190818", "chr6:13190787-13190820",
    "chr6:13190787-13191213", "chr6:19792888-19792974", "chr6:19792891-19792979",
    "chr8:134070215-134070286", "chr8:134070215-134070682", "chr9:112014267-112014291",
    "chr9:112014271-112014290", "chr9:112014271-112014368", "chr9:112014271-112014741",
    "chr9:112014273-112014291", "chrX:11713284-11713308", "chrX:11713284-11713311",
    "chrX:11713284-11713327", "chrX:11713287-11713311", "chrX:11713290-11713318",
    "chrX:11934868-11935072"
  ),
  avrg_meth_promoter = c(
    21.68400000, 0.32142857, 7.48263514, 7.48263514, 7.48263514, 7.48263514,
    3.71646154, 6.78709677, 6.78709677, 6.78709677, 6.78709677, 6.78709677,
    6.78709677, 12.34407895, 0.08561644, 8.34691358, 6.22371795, 6.22371795,
    6.22371795, 6.22371795, 6.22371795, 6.22371795, 6.22371795, 6.22371795,
    6.22371795, 5.72985714, 5.72985714, 6.08114286, 6.08114286, 5.90408451,
    5.90408451, 5.90408451, 5.90408451, 5.90408451, 55.76423913, 55.76423913,
    55.76423913, 55.76423913, 55.76423913, 25.78285714
  )
)

source2_meth_body <- tibble::tibble(
  source2 = c(
    "chr11:24334007-24336603", "chr11:59490054-59491029", "chr2:143253214-143253651",
    "chr2:143253214-143253674", "chr2:143253214-143253679", "chr2:143253215-143253655",
    "chr2:190613314-190613981", "chr22:28669325-28669371", "chr22:28669325-28669921",
    "chr22:28669325-28670131", "chr22:28669325-28670432", "chr22:28669541-28670131",
    "chr22:28669551-28670131", "chr3:123871865-123872261", "chr5:108250797-108251136",
    "chr5:110143384-110144564", "chr6:13190784-13190810", "chr6:13190784-13190820",
    "chr6:13190784-13191209", "chr6:13190784-13191211", "chr6:13190784-13191212",
    "chr6:13190787-13190817", "chr6:13190787-13190818", "chr6:13190787-13190820",
    "chr6:13190787-13191213", "chr6:19792888-19792974", "chr6:19792891-19792979",
    "chr8:134070215-134070286", "chr8:134070215-134070682", "chr9:112014267-112014291",
    "chr9:112014271-112014290", "chr9:112014271-112014368", "chr9:112014271-112014741",
    "chr9:112014273-112014291", "chrX:11713284-11713308", "chrX:11713284-11713311",
    "chrX:11713284-11713327", "chrX:11713287-11713311", "chrX:11713290-11713318",
    "chrX:11934868-11935072"
  ),
  avrg_meth_body = c(
    74.68617, 38.72550, 27.87675, 27.87675, 27.87675, 27.87675, 43.57156,
    76.23071, 76.23071, 76.23071, 76.23071, 76.23071, 76.23071, 44.62959,
    45.98376, 23.61474, 26.59777, 26.59777, 26.59777, 26.59777, 26.59777,
    26.59777, 26.59777, 26.59777, 26.59777, 45.65284, 45.65284, 24.20870,
    24.20870, 44.45126, 44.45126, 44.45126, 44.45126, 44.45126, 54.25862,
    54.25862, 54.25862, 54.25862, 54.25862, 43.26842
  )
)

# 13. Merge methylation data with L1 partners
avg_meth_B42 <- B42_partnered_L1s %>%
  left_join(source2_meth_promoter, by = "source2") %>%
  left_join(source2_meth_body, by = "source2")

# 14. Load and prepare CNV data
B42_CNV <- read_tsv("B42tumor.cov.gz") %>%
  select(chr, start, end, B42tumor_CN) %>%
  dplyr::rename(value = B42tumor_CN) %>%
  mutate(value = log2(value + 1e-3))  # avoid log2(0)

cnv_track <- B42_CNV %>%
  select(chr, start, end, value)


# 15. Initialize Circos plot

# Initialize circos plot
circos.clear()

# Initialize with human ideogram
circos.initializeWithIdeogram(species = "hg38", ideogram.height=convert_height(1, "mm"))

# Add CNV as a line track
circos.genomicTrackPlotRegion(
  data = cnv_track,
  numeric.column = 4,
  track.height = 0.16,
  ylim = c(-2, 2.5),  
  bg.border = "black",
  panel.fun = function(region, value, ...) {
    circos.genomicPoints(
      region = region,
      value = value,
      pch = 16,
      cex = 0.1,
      col = "black"
    )
  }
)

# Prepare methylation heatmap data
heatmap_data_promoter <- avg_meth_B42 %>%
  filter(!is.na(avrg_meth_promoter)) %>%
  mutate(
    chr = gsub(":.*", "", source2),
    start = as.numeric(gsub(".*:(\\d+)-.*", "\\1", source2)),
    end = as.numeric(gsub(".*-(\\d+)", "\\1", source2))
  ) %>%
  select(chr, start, end, avrg_meth_promoter)

# Group by avrg_meth and compute the min and max of start and end for each unique avrg_meth
new_heatmap_data_promoter <- heatmap_data_promoter %>%
  group_by(avrg_meth_promoter) %>%
  summarise(
    chr = chr,  
    start = max(start),
    end = max(end),
    avrg_meth_promoter = unique(avrg_meth_promoter)  # Keep unique avrg_meth value
  ) %>%
  ungroup() %>%
  mutate(
    start = as.numeric(start),
    end = as.numeric(end)
  )%>%
  select(chr,start, end,avrg_meth_promoter) %>%
  arrange(chr, start) %>%
  unique()

# Add heatmap track (colored by methylation value 0–100)
col_fun <- colorRamp2(c(0, 50, 100), c("blue", "white", "red"))

circos.genomicHeatmap(
  new_heatmap_data_promoter,
  heatmap_height = 0.03,
  col = col_fun,
  border = NA,
  side = "inside",connection_height=0.05
)

# Prepare methylation heatmap data
heatmap_data_body <- avg_meth_B42 %>%
  filter(!is.na(avrg_meth_body)) %>%
  mutate(
    chr = gsub(":.*", "", source2),
    start = as.numeric(gsub(".*:(\\d+)-.*", "\\1", source2)),
    end = as.numeric(gsub(".*-(\\d+)", "\\1", source2))
  ) %>%
  select(chr, start, end, avrg_meth_body)

# Group by avrg_meth and compute the min and max of start and end for each unique avrg_meth
new_heatmap_data_body <- heatmap_data_body %>%
  group_by(avrg_meth_body) %>%
  summarise(
    chr = chr,  
    start = max(start),
    end = max(end),
    avrg_meth_body = unique(avrg_meth_body) 
  ) %>%
  ungroup() %>%
  mutate(
    start = as.numeric(start),
    end = as.numeric(end)
  )%>%
  select(chr,start, end,avrg_meth_body) %>%
  arrange(chr, start) %>%
  unique()


# Add heatmap track (colored by methylation value 0–100)
col_fun <- colorRamp2(c(0, 50, 100), c("blue", "white", "red"))

circos.genomicHeatmap(
  new_heatmap_data_body,
  heatmap_height = 0.03,
  col = col_fun,
  border = NA,
  side = "inside",
  connection_height=NULL,
)

B42_partnered_L1s_meth <- avg_meth_B42 %>%
  na.omit()

circos.par(track.margin = c(0.05, 0.05))  # Slight spacing between tracks

# Add links between source and target
for (i in 1:nrow(B42_partnered_L1s_meth)) {
  target_chr <- B42_partnered_L1s_meth$chr[i]
  target_pos <- B42_partnered_L1s_meth$pos[i]
  
  source_info <- strsplit(B42_partnered_L1s_meth$source2[i], "[:-]")[[1]]
  source_chr <- source_info[1]
  source_start <- as.numeric(source_info[2])
  source_end <- as.numeric(source_info[3])
  source_mid <- (source_start + source_end) / 2
  
  # Highlight L1 muultijumps
  is_red <- (
    (source_chr == "chrX" & target_chr == "chr5" & abs(source_mid - 11713304) < 50 & abs(target_pos - 108251482) < 100) |
      (source_chr == "chr5" & target_chr == "chr18" & abs(source_mid - 108250966.5) < 50 & abs(target_pos - 26779223) < 100)
  )
  
  is_blue <- (
    source_chr == "chr11" & target_chr == "chr10" &
      abs(source_mid - 59490541.5) < 100 & abs(target_pos - 17492996) < 100
  )
  
  # Assign color
  link_col <- if (is_red) {
    "red"
  } else if (is_blue) {
    "blue"
  } else {
    "#00000033"
  }
  
  
  # Draw the link
  circos.link(
    sector.index1 = source_chr,
    point1 = source_mid,
    sector.index2 = target_chr,
    point2 = target_pos,
    col = link_col,
    lwd = 1,
    arr.width = 0.1,
    directional = 1
  )
}

# Add legend for methylation gradient
lgd <- ComplexHeatmap::Legend(
  at = c(0, 50, 100),
  col_fun = col_fun,
  title = "Methylation (%)",
  direction = "horizontal"
)

# Draw the legend
draw(lgd, x = unit(0.1, "npc"), y = unit(0.0001, "npc"), just = c("left", "bottom"))

