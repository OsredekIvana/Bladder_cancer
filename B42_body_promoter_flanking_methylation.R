library(dplyr)
library(circlize)
library(tidyr)
library(data.table)
library(GenomicRanges)
library(readxl)
library(fuzzyjoin)
library(readr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(rstatix)

#scienceTheme 
txtFontSize=12
axisFontSize=12
axisTtlFontSize=19
lgdTtlFontSize=19
lgdFontSize=12

scienceTheme=theme(panel.grid.major=element_blank(), 
                   panel.grid.minor=element_blank(), 
                   legend.key=element_blank(), 
                   legend.background=element_blank(), 
                   panel.background = element_blank(), 
                   panel.border=element_blank(),
                   strip.background = element_blank(), 
                   axis.line=element_line(size=0.7, color="black"), 
                   axis.text.x=element_text(size=axisFontSize), 
                   axis.text.y=element_text(size=axisFontSize), 
                   axis.title.x=element_text(size=axisTtlFontSize), 
                   axis.title.y=element_text(size=axisTtlFontSize,angle = 90), 
                   legend.title=element_text(size=lgdTtlFontSize, 
                                             face="bold"),
                   legend.text=element_text(size=lgdFontSize),
                   text=element_text(size=txtFontSize), 
                   strip.text.x=element_text(size=axisTtlFontSize))
theme_set(scienceTheme)

#1. Load Promoter Position Info and Calculate Flanking Regions
B42_detailed_L1 <- read_excel("B42_detailed_L1.xlsx")

B42_flanking_regions <- B42_detailed_L1 %>%
  mutate(
    flank_start_fwd = if_else(strand == "fwd", promotor_start - 501, promotor_end + 1),
    flank_end_fwd = if_else(strand == "fwd", promotor_start - 1, promotor_end + 501),
    flank_start_back = if_else(strand == "fwd", body_end + 1, body_end - 501),
    flank_end_back = if_else(strand == "fwd", body_end + 501, body_end - 1)
  ) %>%
  mutate(
    across(c(flank_start_fwd, flank_end_fwd, flank_start_back, flank_end_back), ~ pmax(.x, 1))
  ) %>%
  select(L1, chr, strand, promotor_start, promotor_end, body_end,
         flank_start_fwd, flank_end_fwd, flank_start_back, flank_end_back)

#2. Read and Combine Methylation Assembly Files
path_to_files <- "/path/l1_source_assembly"
files <- list.files(path_to_files, pattern = "_assembly.pileup.bed.gz", full.names = TRUE)

all_data <- list()
for (file in files) {
  data <- read_tsv(file, col_names = FALSE, comment = "#")
  colnames(data) <- c("chrom", "start", "end", "modified_base", "score", "strand", 
                      "start_incl", "end_incl", "color", "Nvalid_cov", "fraction_modified", 
                      "Nmod", "Ncanonical", "Nother_mod", "Ndelete", "Nfail", "Ndiff", "Nnocall")
  data$file_name <- basename(file)
  all_data <- append(all_data, list(data))
}
combined_data <- rbindlist(all_data) %>%
  mutate(chr = str_extract(file_name, "^chr[0-9XY]+"))

#3. Extract Regions for L1 not present in GRCh38 
get_regions <- function(df, region_start, region_end, region_name) {
  df %>%
    filter(L1 %in% c("notInGRCh38", "somatic")) %>%
    select(chr, strand, !!sym(region_start), !!sym(region_end)) %>%
    mutate(strand = ifelse(strand == "fwd", "+", ifelse(strand == "rev", "-", strand))) %>%
    na.omit() %>%
    inner_join(combined_data, by = c("chr", "strand")) %>%
    filter(start <= !!sym(region_end) & end >= !!sym(region_start)) %>%
    mutate(type = region_name) %>%
    select(-all_of(c(region_start, region_end)))
}

notInGRCh38_promoters_meth <- get_regions(B42_detailed_L1, "promotor_start", "promotor_end", "promoter")
notInGRCh38_body_meth <- get_regions(B42_detailed_L1, "body_start", "body_end", "body")

get_flanking_regions <- function(df, start_col, end_col) {
  df %>%
    filter(L1 %in% c("notInGRCh38", "somatic")) %>%
    select(chr, strand, !!sym(start_col), !!sym(end_col)) %>%
    mutate(strand = ifelse(strand == "fwd", "+", ifelse(strand == "rev", "-", strand))) %>%
    na.omit() %>%
    inner_join(combined_data, by = c("chr", "strand")) %>%
    filter(start <= !!sym(end_col) & end >= !!sym(start_col)) %>%
    mutate(type = "flanking") %>%
    select(-all_of(c(start_col, end_col)))
}

notInGRCh38_flanking_forward_meth <- get_flanking_regions(B42_flanking_regions, "flank_start_fwd", "flank_end_fwd")
notInGRCh38_flanking_back_meth <- get_flanking_regions(B42_flanking_regions, "flank_start_back", "flank_end_back")

# Combine all notInGRCh38 methylation regions
notInGRCh38_regions_meth <- rbind(
  notInGRCh38_flanking_back_meth,
  notInGRCh38_flanking_forward_meth,
  notInGRCh38_body_meth,
  notInGRCh38_promoters_meth
) %>%
  select(chr, start, end, fraction_modified, type) %>%
  mutate(fraction_modified = as.numeric(fraction_modified))

#4. Read InGRCh38 Methylation Bed Files
read_meth_bed <- function(filepath, type_label) {
  fread(filepath, header = FALSE) %>%
    select(V1, V2, V3, V11) %>%
    dplyr::rename(chr = V1, start = V2, end = V3, fraction_modified = V11) %>%
    mutate(type = type_label) %>%
    select(chr, start, end, fraction_modified, type)
}

promoter_meth_inGRCh38_modkit <- read_meth_bed(
  "promoter_methylation_inGRCh38.bed", "promoter")
body_meth_inGRCh38_modkit <- read_meth_bed(
  "body_methylation_inGRCh38.bed", "body")
flank_fwd_meth_inGRCh38_modkit <- read_meth_bed(
  "flanking_fwd_methylation_inGRCh38.bed", "flanking")
flank_back_meth_inGRCh38_modkit <- read_meth_bed(
  "flanking_back_methylation_inGRCh38.bed", "flanking")

#5. Combine All Methylation Data
B42_meth_regions <- rbind(
  notInGRCh38_regions_meth,
  flank_back_meth_inGRCh38_modkit,
  flank_fwd_meth_inGRCh38_modkit,
  body_meth_inGRCh38_modkit,
  promoter_meth_inGRCh38_modkit
) %>%
  mutate(fraction_modified = as.numeric(fraction_modified))

#6. Plot
pastel_colors <- c(
  "promoter" = "#B3EBF2",
  "body" = "#FFD1DC",
  "flanking" = "#70DD99"
)

ggplot(B42_meth_regions, aes(x = type, y = fraction_modified, fill = type)) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(
    method = "wilcox.test", label = "p.signif",
    comparisons = list(
      c("promoter", "body"),
      c("promoter", "flanking"),
      c("body", "flanking")
    )
  ) +
  scale_fill_manual(values = pastel_colors) +
  labs(
    title = "Methylation levels",
    x = "Region Type",
    y = "Fraction Modified (%)"
  ) +
  theme(legend.position = "none") +
  scienceTheme

#7. Statistical Tests
pairwise_results <- B42_meth_regions %>%
  wilcox_test(fraction_modified ~ type, p.adjust.method = "BH") %>%
  add_significance()
