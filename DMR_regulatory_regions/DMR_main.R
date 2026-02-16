#Script to analyse dmr of promoter regions between the two haploypes of the same tumour for all samples, and to look for the potential SVs present in the 
#regions and their downstream effects 

library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(data.table)
library(GenomicRanges)
library(stringr)
library(fuzzyjoin)
library(stringr)
library(pheatmap)
library(ggrepel)
library(maftools)



#scienceTheme 
txtFontSize=16
axisFontSize=16
axisTtlFontSize=22
lgdTtlFontSize=22
lgdFontSize=16

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


file_path <- "bed_dmr_files/dmr_regulatory_regions/"

#Read all files that are in this folder 
#Each folder contains corresponding *tumor_hp_DMR file to be read
#Read all files into one  dataframe with the name column containing the corresponding name sample 

# Path to the parent folder containing *tumor_hp subfolders
parent_dir <- "dmr_regulatory_regions/dna_based"

#Find all files ending with "_tumor_hp_DMR" inside subfolders
dmr_files <- list.files(
  path = parent_dir,
  pattern = "tumor_hp_DMR$",
  recursive = TRUE,
  full.names = TRUE)

read_dmr_file <- function(path) {
  
  # Extract sample name from parent folder:
  # e.g., dna_based/B42tumor_hp/B42tumor_hp_DMR → "B42"
  folder <- basename(dirname(path))
  sample <- str_replace(folder, "tumor_hp$", "")
  
  # Read file
  df <- fread(path, header = FALSE)
  
  # Set the correct column names according to modkit DMR output
  colnames(df) <- c(
    "chrom", "start", "end", "name", "score",
    "sampleA_counts", "sampleA_total",
    "sampleB_counts", "sampleB_total",
    "sampleA_percents", "sampleB_percents",
    "sampleA_fraction_modified", "sampleB_fraction_modified"
  )
  
  # Annotate with sample name
  df$sample <- sample
  
  return(df)
}

dmr_all <- rbindlist(lapply(dmr_files, read_dmr_file))

# Preview the result
head(dmr_all)


#Keep only non NA rows 

dmr_all_clean <- dmr_all %>%
  filter(
    !is.na(sampleA_fraction_modified),
    !is.na(sampleB_fraction_modified)
  )

dmr_all_clean <- dplyr::filter(dmr_all_clean, !is.na(score))


#Inspect the scores 
ggplot(dmr_all_clean, aes(x = score)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black") +
  theme_bw() +
  labs(
    title = "Distribution of DMR Score",
    x = "Score",
    y = "Number of promoters"
  )+
  facet_grid(~sample)


#Filter the score 
dmr_all_filtered <- dmr_all_clean %>%
filter(score>20)

ggplot(dmr_all_filtered, aes(x = score)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black") +
  theme_bw() +
  labs(
    title = "Distribution of DMR Score",
    x = "Score",
    y = "Number of promoters"
  )+
  facet_wrap(~sample)



#Compute the effect size
dmr_all_clean$delta_meth <- dmr_all_clean$sampleA_fraction_modified -
  dmr_all_clean$sampleB_fraction_modified

ggplot(dmr_all_clean, aes(x = delta_meth)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black") +
  theme_bw() +
  labs(
    title = "Distribution of Effect Size",
    x = "Effect Size",
    y = "Number of promoters"
  )+
  facet_wrap(~sample)

dmr <- dmr_all_clean[
  is.finite(sampleA_fraction_modified) &
    is.finite(sampleB_fraction_modified) &
    sampleA_total > 0 & sampleB_total > 0
]


#Test whether the difference between modified fractions for each promoter is the significant (perc modified)
dmr_all_with_q <- dmr %>%
  group_by(sample) %>%                  # run per sample
  mutate(
    a_mod = as.numeric(gsub("m:", "", sampleA_counts)),
    b_mod = as.numeric(gsub("m:", "", sampleB_counts)),
    # compute p-value per row using rowwise within mutate
    pvalue = purrr::pmap_dbl(
      list(a_mod, sampleA_total, b_mod, sampleB_total),
      ~ binom.test(c(..1, ..3), c(..2, ..4))$p.value
    ),
    # adjust for multiple testing within this sample
    padj = p.adjust(pvalue, method = "BH")
  ) %>%
  ungroup() %>%
  dplyr::select(-a_mod, -b_mod)

#Final QC filtering
sig_q <- dmr_all_with_q %>%
  filter(abs(delta_meth) > 0.2 &   # 20% diff methylation  
           padj < 0.05 &
           abs(score) > 50) #score is a coverage-weighted difference


# Visualise 
samples <- unique(sig_q$sample)

for(s in samples){
  sig_sample <- sig_q %>% filter(sample == s)
  
  # convert to numeric matrix
  mat_sample <- as.matrix(sig_sample[, c("sampleA_fraction_modified", "sampleB_fraction_modified")])
  
  row_labels <- sig_sample$chrom
  
  pheatmap(mat_sample,
           labels_row = row_labels,
           main = paste0("Differential promoters: ", s))
}




###########
#Intersect all with the variants from the maf file 

#Extract the maf data into a data frame 

summary <- "summary.maf"
meta <- "combined_meta.tsv"
maf = read.maf(maf=summary, clinicalData = meta)
maf_variants <- as.data.frame(maf@data) %>%
  select(Tumor_Sample_Barcode,Hugo_Symbol,Chromosome,Start_Position,End_Position,
         Strand, Variant_Classification)

maf_variants$Tumor_Sample_Barcode <- gsub("tumor$", "", maf_variants$Tumor_Sample_Barcode, ignore.case = TRUE)



#Match the variants 
#Convert each table into genomics ranges 
sig_q <- as.data.frame(sig_q)

hits_dplyr_q <- fuzzy_inner_join(
  sig_q,
  maf_variants,
  by = c(
    "sample" = "Tumor_Sample_Barcode",
    "chrom"  = "Chromosome",
    "start"  = "Start_Position",
    "end"    = "End_Position"
  ),
  # condition for overlap:
  # sig.start <= variant.end  AND  sig.end >= variant.start
  match_fun = list(`==`,`==`, `<=`, `>=`)
)


#Use Genomics Rangens otherwise it is very slow with fizzy join
sig_gr <- GRanges(
  sample = sig_q$sample,
  seqnames = sig_q$chrom,
  ranges   = IRanges(sig_q$start, sig_q$end),
  name     = sig_q$name,
  delta_meth = sig_q$delta_meth,
  padj     = sig_q$padj
)


mut_gr <- GRanges(
  sample = maf_variants$Tumor_Sample_Barcode,
  seqnames = maf_variants$Chromosome,
  ranges   = IRanges(maf_variants$Start_Position,
                     maf_variants$End_Position),
  gene     = maf_variants$Hugo_Symbol,
  variant  = maf_variants$Variant_Classification
)


### Loop over samples
samples <- unique(sig_gr$sample)
hits_list <- lapply(samples, function(s){
  sig_sub  <- sig_gr[sig_gr$sample == s]
  mut_sub  <- mut_gr[mut_gr$sample == s]
  
  # find overlaps
  hits <- findOverlaps(sig_sub, mut_sub)
  
  # return the overlapping rows
  cbind(
    as.data.frame(sig_sub[queryHits(hits)]),
    as.data.frame(mut_sub[subjectHits(hits)])
  )
})

# Combine all hits
hits_all <- do.call(rbind, hits_list)



######## Create volcano plots for B94 with COL4A1 and COL4A2 and B187 and NDUFA1
########  Colored by color gradient for all significant genes that reflects the normalized gene expression #######

dmr_B94  <- dmr_all_with_q %>% filter(sample == "B94")


#1. dmr regions already have a column name that corresponds to NCBI gene id. 
library(AnnotationDbi)
library(org.Hs.eg.db)


entrez_to_ensembl <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys     = as.character(unique(dmr_B94$name)),
  keytype  = "ENTREZID",
  columns  = c("ENSEMBL", "SYMBOL")
)

dmr_B94 <- dmr_B94 %>%
  mutate(name = as.character(name))  %>%
  left_join(
    entrez_to_ensembl,
    by = c("name" = "ENTREZID")
  )

#2. Read in the normalized gene expression from salmon 

expr_B94 <- read.delim("dmr_regulatory_regions/quant.genes94.sf")
expr_B94 <- expr_B94 %>%
  dplyr::select(
    ENSEMBL = Name,
    TPM
  ) %>%
  mutate(
    logTPM = log2(TPM + 1)
  )

dmr_B94 <- dmr_B94 %>%
  left_join(expr_B94, by = "ENSEMBL")

#Sanity check 
# How many DMRs have expression info?
dmr_B94 %>%
  summarise(
    total = n(),
    with_expression = sum(!is.na(TPM))
  )

# Check genes of interest
dmr_B94 %>%
  filter(SYMBOL %in% c("COL4A1", "COL4A2")) %>%
  dplyr::select(SYMBOL, ENSEMBL, TPM, logTPM) %>%
  distinct()

#Filter 
dmr_B94_gene <- dmr_B94 %>%
  mutate(
    significant =
      abs(delta_meth) > 0.2 &
      padj < 0.05 &
      score > 50
  )

#Sanity check

table(dmr_B94_gene$significant)

dmr_B94_gene %>%
  filter(SYMBOL %in% c("COL4A1", "COL4A2")) %>%
  dplyr::select(SYMBOL, delta_meth, padj, score, logTPM, significant)


#Add COSMIC genes  (list obtained from https://cancer.sanger.ac.uk/cosmic/download/cosmic/v103/cancergenecensus)

#Add to the main table 
dmr_B94_gene <- dmr_B94_gene %>%
  mutate(is_cosmic = SYMBOL %in% cosmic_gene_list)

highlight_cosmic_sig <- dmr_B94_gene %>%
  filter(is_cosmic & significant)


#Volcano plot visualization 
# Genes to highlight
highlight_genes <- dmr_B94_gene %>% filter(SYMBOL %in% c("COL4A2"))

ggplot(dmr_B94_gene, aes(
  x = delta_meth,
  y = -log10(padj)
)) +
  # non-significant genes
  geom_point(
    data = dmr_B94_gene %>% filter(!significant),
    color = "grey80",
    size = 1
  ) +
  
  # significant genes, coloured by expression
  geom_point(
    data = dmr_B94_gene %>% filter(significant),
    aes(color = logTPM),
    size = 1.5
  ) +
  
  scale_color_viridis_c(
    option = "plasma",
    name = "log2(TPM + 1)",
    na.value = "grey80"
  ) +
  
  # Circle the highlighted gene(s)
  geom_point(
    data = highlight_genes,
    aes(x = delta_meth, y = -log10(padj)),
    shape = 21,       # hollow circle
    color = "red",
    size = 3,
    stroke = 1      # thickness of the circle border
  ) +
  
  # Label the gene(s)
  geom_text_repel(
    data = highlight_genes,
    aes(label = SYMBOL),
    nudge_x = 0.05,
    nudge_y = 1,
    size = 4,
    box.padding = 0.3
  ) +
  
  # Threshold lines
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  
  theme_bw() +
  labs(
    x = "Mean haplotype methylation difference",
    y = "-log10(adjusted p-value)"
  ) +
  
  geom_point(
    data = highlight_cosmic_sig,
    shape = 21, color = "black", size = 3, stroke = 1
  ) +
  geom_text_repel(
    data = highlight_cosmic_sig,
    aes(label = SYMBOL), size = 4
  ) +
  scienceTheme


########################## Repeat for the B187 ##########################



#Repeat for the B187
dmr_B187  <- dmr_all_with_q %>% filter(sample == "B187")


entrez_to_ensembl <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys     = as.character(unique(dmr_B187$name)),
  keytype  = "ENTREZID",
  columns  = c("ENSEMBL", "SYMBOL")
)

dmr_B187 <- dmr_B187 %>%
  mutate(name = as.character(name))  %>%
  left_join(
    entrez_to_ensembl,
    by = c("name" = "ENTREZID")
  )

#2. Read in the normalized gene expression from salmon 

expr_B187 <- read.delim("dmr_regulatory_regions/quant.genes187.sf")
expr_B187 <- expr_B187 %>%
  dplyr::select(
    ENSEMBL = Name,
    TPM
  ) %>%
  mutate(
    logTPM = log2(TPM + 1)
  )

dmr_B187 <- dmr_B187 %>%
  left_join(expr_B187, by = "ENSEMBL")

#Sanity check 
# How many DMRs have expression info?
dmr_B187 %>%
  summarise(
    total = n(),
    with_expression = sum(!is.na(TPM))
  )

# Check genes of interest
dmr_B187 %>%
  filter(SYMBOL %in% c("NDUFA1")) %>%
  dplyr::select(SYMBOL, ENSEMBL, TPM, logTPM) %>%
  distinct()

#Filter 
dmr_B187_gene <- dmr_B187 %>%
  mutate(
    significant =
      abs(delta_meth) > 0.2 &
      padj < 0.05 &
      score > 50
  )

#Sanity check

table(dmr_B187_gene$significant)

dmr_B187_gene %>%
  filter(SYMBOL %in% c("NDUFA1")) %>%
  dplyr::select(SYMBOL, delta_meth, padj, score, logTPM, significant)


#Add cosmic genes 


#Add to the main table 
dmr_B187_gene <- dmr_B187_gene %>%
  mutate(is_cosmic = SYMBOL %in% cosmic_gene_list)

highlight_cosmic_sig <- dmr_B187_gene %>%
  filter(is_cosmic & significant)


#Volcano plot visualization 
# Genes to highlight
highlight_genes <- dmr_B187_gene %>% filter(SYMBOL %in% c("NDUFA1"))

ggplot(dmr_B187_gene, aes(
  x = delta_meth,
  y = -log10(padj)
)) +
  # non-significant genes
  geom_point(
    data = dmr_B187_gene %>% filter(!significant),
    color = "grey80",
    size = 1
  ) +
  
  # significant genes, coloured by expression
  geom_point(
    data = dmr_B187_gene %>% filter(significant),
    aes(color = logTPM),
    size = 1.5
  ) +
  
  scale_color_viridis_c(
    option = "plasma",
    name = "log2(TPM + 1)",
    na.value = "grey80"
  ) +
  
  # Circle the highlighted gene(s)
  geom_point(
    data = highlight_genes,
    aes(x = delta_meth, y = -log10(padj)),
    shape = 21,       # hollow circle
    color = "red",
    size = 3,
    stroke = 1      # thickness of the circle border
  ) +
  
  # Label the gene(s)
  geom_text_repel(
    data = highlight_genes,
    aes(label = SYMBOL),
    nudge_x = 0.05,
    nudge_y = 1,
    size = 4,
    box.padding = 0.3
  ) +
  
  # Threshold lines
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  
  theme_bw() +
  labs(
    x = "Mean haplotype methylation difference",
    y = "-log10(adjusted p-value)"
  ) +
  
  geom_point(
    data = highlight_cosmic_sig,
    shape = 21, color = "black", size = 3, stroke = 1
  ) +
  geom_text_repel(
    data = highlight_cosmic_sig,
    aes(label = SYMBOL), size = 4
  ) +
  scienceTheme


