library(dplyr)
library(ComplexHeatmap)

literature_drivers <- c("TP53", "FGFR3", "CREBBP", "EP300","ERBB3","ERCC2","HRAS","KDM6A","KMT2D","PIK3CA",
                        "KMT2A", "SPTAN1", "ERBB2", "CREBBP", "FAT1", "ATM", "KMT2C",
                        "P53", "PIK3CA", "TSC1", "FGFR3", "HRAS", "HER2",
                        "P16","EGFR", "ERBB2", "ERBB3","SPATC1L","NAPRT","PARP6",
                        "LXN","RB1", "NF2", "NF1", "TSC2", "TSC1", "PTEN",
                        "ERCC1", "XRCC1", "GSTP1", "CDA", "GSTM1","GSTT","TERT",
                        "UTX", "MLL-MLL3", "CREBBP-EP300", "NCOR1", "ARID1A", "CHD6",
                        "PD-L1","BRCA1", "STAG2","ESPL1","KANSL1","ASXL2","LAPTM4B",
                        "ARHGAP35","KLF5","FBXW7","CDKN2A","PTEN","EP300","SPTAN1",
                        "CDKN1A","RAI1","RHOA","FGFR3","NRAS","ELF3","HRAS",
                        "OTOP1","ZFP36L1","FRG1","UQCR10","RHOB","HES1","PIGP")

#Calculate the percentage of samples in which a given gene is mutated
total_samples <- cn.data %>% 
  pull(Sample_name) %>% 
  n_distinct()

# for each gene, count how many unique samples have Amp or Del, then turn into a percentage
gene_perc <- cn.data %>%
  filter(CN %in% c("Amp", "Del")) %>%
  group_by(Gene) %>%
  summarise(
    n_aberrant = n_distinct(Sample_name),
    pct_aberrant = n_aberrant / total_samples * 100
  ) %>%
  arrange(desc(pct_aberrant))

gene_perc

#filter gene_perc for the literature_drivers
driver_perc <- gene_perc %>%
  filter(Gene %in% literature_drivers)


#filter the original data for the  literature_drivers that are modified in more than 8% of samples  (based on the previous table)

# Identify literature‐driver genes with > 8% aberration
high_freq_drivers <- driver_perc %>%
  filter(pct_aberrant > 8) %>%
  pull(Gene)

# Filter the original CNV calls for those high-frequency drivers
filtered_cn_data <- cn.data %>%
  filter(
    Gene %in% high_freq_drivers)


# Define your sample list and driver genes 
samples <- c("B4tumor","B060tumor","B175tumor","B42tumor","B134tumor","B5tumor",
             "B157tumor","B67tumor","B123tumor","B187tumor","B74tumor","B32tumor",
             "B85tumor","B12tumor","B156tumor","B125tumor","B24tumor","B154tumor",
             "B39tumor","B22tumor","B87tumor","B94tumor","B178tumor")


# Build a binary presence/absence matrix: rows = genes, columns = samples
mat <- matrix(0,
              nrow = length(high_freq_drivers),
              ncol = length(samples),
              dimnames = list(high_freq_drivers, samples))

#Fill in 1 where cn.data indicates Amp or Del
filtered_cn_data %>%
  filter(Gene %in% high_freq_drivers, Sample_name %in% samples) %>%
  distinct(Gene, Sample_name) %>%
  rowwise() %>%
  do({
    mat[.$Gene, .$Sample_name] <<- 1
    data.frame()  # dummy
  })

