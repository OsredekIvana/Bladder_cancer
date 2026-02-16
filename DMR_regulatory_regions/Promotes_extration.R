BiocManager::install("TxDb.Hsapiens.UCSC.hg38.refGene")

library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(GenomicFeatures)
library(GenomicRanges)
library(dplyr)


# Load the TxDb object
txdb <- TxDb.Hsapiens.UCSC.hg38.refGene

# Extract promoters: 2000 bp upstream + 400 bp downstream
promoters_gr <- promoters(txdb, upstream = 2000, downstream = 200)

# Convert to a plain dataframe without row names
promoters_df_transcript <- as.data.frame(promoters_gr, row.names = NULL)

# Remove exact duplicate coordinates
promoters_unique <- promoters_df %>%
  distinct(seqnames, start, end, strand, .keep_all = TRUE)

promoters_filtered <- promoters_unique %>%
  filter(!grepl("alt", seqnames)) %>%
  filter(!grepl("fix", seqnames)) %>%
  filter(!grepl("chrUn", seqnames)) %>%
  filter(!grepl("random", seqnames))

bed <- promoters_filtered %>%
  dplyr::select(seqnames, start, end, tx_name)

write.table(bed, "mdr_regulatory_regions/hg38_promoters.bed",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)



###### based on genes rather than transcripts 
genes_gr <- genes(txdb)

# For one promoter per gene, use promoters() on genes
upstream <- 2000
downstream <- 400
proms_by_gene <- promoters(genes_gr, upstream=upstream, downstream=downstream)
promoters_df_dna <- as.data.frame(proms_by_gene, row.names = NULL)

# Remove exact duplicate coordinates
promoters_unique_dna <- promoters_df_dna %>%
  distinct(seqnames, start, end, strand, .keep_all = TRUE)

promoters_filtered_dna <- promoters_unique_dna %>%
  filter(!grepl("alt", seqnames)) %>%
  filter(!grepl("fix", seqnames)) %>%
  filter(!grepl("chrUn", seqnames)) %>%
  filter(!grepl("random", seqnames))


bed <- promoters_filtered_dna %>%
  dplyr::select(seqnames, start, end, gene_id)

chrom_order <- paste0("chr", c(1:22, "X", "Y", "M"))

bed_sorted <- bed %>%
  mutate(seqnames = factor(seqnames, levels = chrom_order)) %>%
  arrange(seqnames, start)



write.table(bed_sorted, "dmr_regulatory_regions/hg38_promoters_dna.bed",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
