library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)

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


#1. Read the L1 data 
L1_data <- read_excel("suppl_table3.xlsx", col_names = TRUE)  

#2. Select relevant columns 
L1_filtered <- L1_data %>%
  select(tumor_wgs_icgc_sample_id,histology_abbreviation,nb_L1,nb_L1_solo,nb_TD,nb_DEL,nb_DUP)

#3. Read the data containing tumour_aliqote ID in pcawg 
pcawg_info <- read.table("release_may2016.v1.2.tsv", header=T, sep="\t")
pcawg_info <- select(pcawg_info, donor_unique_id,dcc_project_code,tumor_wgs_aliquot_id,tumor_wgs_icgc_sample_id)

#4. Filter rows with multiple entries (comma separate entries in the same row of tumor_wgs_icgc_sample_id columns)
pcawg_info_filtered <- pcawg_info %>%
  filter(grepl(",", tumor_wgs_icgc_sample_id))

#5. Remove those rows from the original pcawg_info dataset
pcawg_info_cleaned <- pcawg_info %>%
  filter(!grepl(",", tumor_wgs_icgc_sample_id))

#6. Amplicon/CNV metadata 
metadata <- read.table("tumour_cnv_data_p53.tsv", header=T)

#7. Merge the  tables based on the tumour ids
pcawg_L1 <- full_join(L1_filtered,pcawg_info_cleaned, join_by(tumor_wgs_icgc_sample_id == tumor_wgs_icgc_sample_id))

#8. Remove rows with NA in donor_unique_id
pcawg_L1_cleaned <- pcawg_L1 %>%
  filter(!is.na(donor_unique_id))

L1_metadata <- left_join(pcawg_L1_cleaned,metadata,join_by(donor_unique_id == donor_unique_id))
L1_metadata <- L1_metadata %>%
  mutate(amplicon = ifelse(CNV_type=="AMPLICON", "yes", "no")) %>%
  na.omit()

#9. Compute p-value
pval <- compare_means(nb_L1 ~ TP53_mutations, data = L1_metadata, method = "wilcox.test")$p.format

#10. Plot 
ggplot(L1_metadata, aes(x = amplicon, y = nb_L1, fill = TP53_mutations)) +
  geom_boxplot() +                                       
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red")) + 
  scale_y_log10() +
  stat_compare_means(method = "wilcox.test", label = "p.signif") +
  labs(
    subtitle = paste("Wilcoxon test p-value:", pval),
    fill = paste("TP53 Mutation")
  )





