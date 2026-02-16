library(GenomicRanges)
library(BSgenome)
library(MutationalPatterns)
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(patchwork)


#Recreating the figure with clonal/subclonal ratios for the paper using the standard theme 
txtFontSize=15
axisFontSize=15
axisTtlFontSize=20
lgdTtlFontSize=20
lgdFontSize=15

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

# Function to create GenomicRanges object
GRmaker = function(CALLS, GROUP) {
  TEMP = CALLS[CALLS$GROUP == GROUP,]
  GR = with(TEMP, GRanges(CHROM, ranges=IRanges(POS, POS), REF=REF, ALT=ALT, VAF=T_VAF))
  seqlevelsStyle(GR) = "UCSC"
  GenomeInfoDb::genome(GR) = "hg38"
  return(GR)
}

# Load the data
df = read.table("/Users/osredek/Desktop/EMBL/PROJECT_Bladder_cancer/Clonality_signatures/Revisions/mutations.tsv", header=T, sep='\t')

clonal = df %>% filter(T_VAF > 0.3)
subclonal = df %>% filter(T_VAF < 0.3)

# Load reference genome
ref_genome = "BSgenome.Hsapiens.UCSC.hg38"
colors = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928", "#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F", "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99")


# Add GROUP column to clonal and subclonal data
clonal$GROUP = clonal$sampleID
subclonal$GROUP = subclonal$sampleID

# Get unique sample groups
GROUPS_clonal = unique(clonal$GROUP)
GROUPS_subclonal = unique(subclonal$GROUP)

# Create GRmaker lists for each group
G_lists_clonal = lapply(GROUPS_clonal, function(x) GRmaker(clonal, x))
names(G_lists_clonal) = GROUPS_clonal

G_lists_subclonal = lapply(GROUPS_subclonal, function(x) GRmaker(subclonal, x))
names(G_lists_subclonal) = GROUPS_subclonal  

# Calculate mutation type occurrences
type_occurences_clonal = mut_type_occurrences(G_lists_clonal, ref_genome)
type_occurences_subclonal = mut_type_occurrences(G_lists_subclonal, ref_genome)

# Plot the SNV spectrum and save as PNG
p = plot_spectrum(type_occurences_clonal, indv_points=T, by=rownames(type_occurences_clonal), error_bars="none")
p = p + scienceTheme
p
dev.off()

# Create mutation matrices
mut_mat_clonal = mut_matrix(vcf_list = G_lists_clonal, ref_genome = ref_genome)
png("/Users/osredek/Desktop/EMBL/PROJECT_Bladder_cancer/Clonality_signatures/Revisions/snv_context_clonal.png", width=1200, height=600)
p=plot_96_profile(mut_mat_clonal)
p=p + scienceTheme
dev.off()  # Close the PNG device


plot_spectrum(type_occurences_subclonal, indv_points=T, by=rownames(type_occurences_subclonal), error_bars="none")
mut_mat_subclonal = mut_matrix(vcf_list = G_lists_subclonal, ref_genome = ref_genome)
png("/Users/osredek/Desktop/EMBL/PROJECT_Bladder_cancer/Clonality_signatures/Revisions/snv_context_subclonal.png", width=1200, height=600)
plot_96_profile(mut_mat_subclonal)
dev.off()  # Close the PNG device


# Load mutational signatures
signatures = get_known_signatures()
signatures = signatures[,c("SBS1","SBS2","SBS3","SBS4","SBS5","SBS7b","SBS8","SBS13","SBS25","SBS26","SBS30","SBS37","SBS41","SBS92","SBS93","SBS94")]

# Fit mutational signatures to the data
refit_clonal = fit_to_signatures_strict(mut_mat_clonal, signatures)
refit_subclonal = fit_to_signatures_strict(mut_mat_subclonal, signatures)

# Extract and save contribution tables
ctrb_clonal = as.data.frame(refit_clonal$fit_res$contribution)
ctrb_clonal = sweep(ctrb_clonal, 2, colSums(ctrb_clonal), `/`)
ctrb_clonal$signature = rownames(ctrb_clonal)
write.table(ctrb_clonal, file="/Users/osredek/Desktop/EMBL/PROJECT_Bladder_cancer/Clonality_signatures/Revisions/clonal_signatures.tsv", sep="\t", quote=FALSE, row.names=FALSE)

ctrb_subclonal = as.data.frame(refit_subclonal$fit_res$contribution)
ctrb_subclonal = sweep(ctrb_subclonal, 2, colSums(ctrb_subclonal), `/`)
ctrb_subclonal$signature = rownames(ctrb_subclonal)
write.table(ctrb_subclonal, file="/Users/osredek/Desktop/EMBL/PROJECT_Bladder_cancer/Clonality_signatures/Revisions/subclonal_signatures.tsv", sep="\t", quote=FALSE, row.names=FALSE)


# Plot and save contributions
png("/Users/osredek/Desktop/EMBL/PROJECT_Bladder_cancer/Clonality_signatures/Revisions/clonal_signatures_rel.png", width=1200, height=600)
plot_contribution(refit_clonal$fit_res$contribution, coord_flip=TRUE, mode="relative")
dev.off()  # Close the PNG device

png("/Users/osredek/Desktop/EMBL/PROJECT_Bladder_cancer/Clonality_signatures/Revisions/clonal_signatures_abs.png", width=1200, height=600)
plot_contribution(refit_clonal$fit_res$contribution, coord_flip=TRUE, mode="absolute")
dev.off()  # Close the PNG device

plot_contribution(refit_subclonal$fit_res$contribution, coord_flip=TRUE, mode="relative")
png("/Users/osredek/Desktop/EMBL/PROJECT_Bladder_cancer/Clonality_signatures/Revisions/subclonal_signatures_rel.png", width=1200, height=600)
dev.off()  # Close the PNG device

plot_contribution(refit_subclonal$fit_res$contribution, coord_flip=TRUE, mode="absolute")
png("/Users/osredek/Desktop/EMBL/PROJECT_Bladder_cancer/Clonality_signatures/Revisions/subclonal_signatures_abs.png", width=1200, height=600)
dev.off()  # Close the PNG device

# Clonal plot (with labels and legend)
p1 <- plot_contribution(refit_clonal$fit_res$contribution, coord_flip = TRUE, mode = "relative") + 
  ggtitle("Clonal") + scienceTheme


# Subclonal plot (no sample labels, no legend)
p2 <- plot_contribution(refit_subclonal$fit_res$contribution, coord_flip = TRUE, mode = "relative") + 
  ggtitle("Subclonal") + scienceTheme

# Combine them tightly with shared title
(p1 + p2) + 
  plot_layout(widths = c(1, 1), guides = "collect")   



# Plot and save reconstruction error
png("/Users/osredek/Desktop/EMBL/PROJECT_Bladder_cancer/Clonality_signatures/Revisions/clonal_reconstruction_error.png", width=1200, height=600)
plot_original_vs_reconstructed(mut_mat_clonal, refit_clonal$fit_res$reconstructed, y_intercept=0.9)
dev.off()  # Close the PNG device

png("/Users/osredek/Desktop/EMBL/PROJECT_Bladder_cancer/Clonality_signatures/Revisions/subclonal_reconstruction_error.png", width=1200, height=600)
plot_original_vs_reconstructed(mut_mat_subclonal, refit_subclonal$fit_res$reconstructed, y_intercept=0.9)
dev.off()  # Close the PNG device


##Calculate ratio of clonal vs subclonal signatures for each tumour: x-axis plot the signatures, y-axis the ratio and the value as the boxplot for each tumor
#1. Turn the contributions to relative contribution (0-1)
#2. Calculate the ratio between clonal vs subclonal
#3. x-axsis plot the signatures using boxplot, y-axis ratio 

ctrb_subclonal
ctrb_clonal


#Subset ctrb_subclonal and ctrb_clonal only to tumours with a cosine similarity > 0.95 


#B74 is missing from the original SNV data!!
ctrb_subclonal_QCed <- ctrb_subclonal %>%
  select("signature","B060","B123","B125","B12","B134",
         "B156","B157","B175","B178","B187","B32","B24",
         "B39",
         "B42","B4", "B5", "B67","B85",
         "B94") %>%
  tibble::remove_rownames()

#excluding B22 because it has no subclonal mutations 
ctrb_clonal_QCed <- ctrb_clonal %>%
  select("signature","B060","B123","B125","B12","B134",
         "B32",
         "B156","B157","B175","B178","B187","B24",
         "B39",
         "B42","B4", "B5", "B67","B85",
         "B94")%>%
  tibble::remove_rownames()

#1. Calculate relative contribution of each signature in a particular tumour 

ctrb_clonal_QCed_relative <-ctrb_clonal_QCed

for (col in colnames(ctrb_clonal_QCed_relative)[-1]) {
  total <- sum(ctrb_clonal_QCed_relative[[col]])  # Calculate the column sum
  ctrb_clonal_QCed_relative[[col]] <- (ctrb_clonal_QCed_relative[[col]] / total)  
}


ctrb_subclonal_QCed_relative <-ctrb_subclonal_QCed

for (col in colnames(ctrb_subclonal_QCed_relative)[-1]) {
  total <- sum(ctrb_subclonal_QCed_relative[[col]])  # Calculate the column sum
  ctrb_subclonal_QCed_relative[[col]] <- (ctrb_subclonal_QCed_relative[[col]] / total)  
}

#2. Calculate the ratio between clonal vs subclonal

ctrb_subclonal_QCed_relative <-ctrb_subclonal_QCed_relative %>%
  mutate(status="subclonal")

ctrb_clonal_QCed_relative <-ctrb_clonal_QCed_relative %>%
  mutate(status="clonal")

# Initialize an empty dataframe to store the ratio
ctrb_ratio_QCed <- ctrb_clonal_QCed_relative

# Iterate over all tumor columns (excluding "signature" and "status")
for (col in colnames(ctrb_ratio_QCed)[-c(1, ncol(ctrb_ratio_QCed))]) {
  ctrb_ratio_QCed[[col]] <- ctrb_clonal_QCed_relative[[col]] / ctrb_subclonal_QCed_relative[[col]]
}
ctrb_ratio_QCed$status <- "clonal/subclonal ratio"

#3. plot: x-axsis  the signatures using boxplot, y-axis ratio 
ctrb_ratio_long <- ctrb_ratio_QCed %>%
  pivot_longer(cols = starts_with("B"), 
               names_to = "tumor", 
               values_to = "ratio")

# Create the boxplot
ggplot(ctrb_ratio_long, aes(x = signature, y = ratio)) +
  geom_boxplot() +
  #theme_minimal() +
  geom_hline(yintercept = 1, color = "red") + 
  labs(x = "Signature", y = "Clonal/Subclonal Ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  
  

#Adding smoking info to the signatures 
# Smoking status data
smoking_status <- tibble(
  Sample = c("B123", "B125", "B12", "B134", "B157", "B178", "B24", "B32", "B39", "B42", 
             "B4", "B67", "B74", "B85", "B94", "B60", "B154", "B156", "B175", "B187", 
             "B22", "B5", "B87"),
  smoking_status = c("y", "y", "y", "y", "n", "y", "y", "y", "y", "y", 
                     "y", "y", "y", "y", "y", "n", "n", "n", "n", "n", 
                     "n", "n", "n")
)

# Extract sample ID from tumor column
ctrb_ratio_long <- ctrb_ratio_long %>%
  mutate(Sample = gsub("tumor", "", tumor))

# Merge the tables
ctrb_ratio_long <- ctrb_ratio_long %>%
  left_join(smoking_status, by = "Sample")

SBS92 <- ctrb_ratio_long %>%
  filter(signature == "SBS92")

ggplot(SBS92, aes(x = signature, y = ratio, fill = signature)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Signature", y = "Clonal/Subclonal Ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  + facet_wrap(~smoking_status)




#Remove signatures for which the effect is unknown: SBS25, SBS37, SBS41, SBS94 

ctrb_ratio_long_filtered <- ctrb_ratio_long %>%
  filter(!signature %in% c("SBS25", "SBS37", "SBS41", "SBS94","SBS7b","SBS8","SBS93")) %>%
  #remove SBS26 that has no datapoints 
  filter(!signature == "SBS26")


# Create the boxplot
ggplot(ctrb_ratio_long_filtered, aes(x = signature, y = ratio)) +
  stat_boxplot(geom = "errorbar")+
  geom_boxplot() +
  geom_hline(yintercept = 1, color = "red") + 
  labs(x = "Signature", y = "Clonal/Subclonal Ratio") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.3, vjust = 0.5),
        axis.title.y = element_text(vjust = 2))

