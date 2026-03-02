################################################################################
# Placental Long-Read RNA-seq Transcriptome Assembly Analysis
# Author: Sean T. Bresnahan
# Description: Analysis of placental transcriptome assembly against GENCODEv45
#              reference using SQANTI3 classification and comparison with GTEx v9
#
# Input Datasets:
#   - Supplementary/ONT_RIN_scores.csv
#       RNA Integrity Number (RIN) scores for Oxford Nanopore (ONT) long-read samples
#   - Supplementary/GUSTO_RIN_scores.csv
#       RIN scores for GUSTO short-read validation cohort
#   - Assembly/ESPRESSO_corrected_classification.txt
#       SQANTI3 structural classification of assembled transcripts
#   - Assembly/ESPRESSO_corrected_corrected.gtf
#       GTF annotation file for corrected long-read assembly
#   - Assembly/ESPRESSO_corrected_SC_filtered_classification.txt
#       SQANTI3 classification after structural category filtering
#   - /rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/GenomicReferences/txome/gencode_v45/gencode.v45.annotation.gtf
#       GENCODE v45 reference annotation
#   - /rsrch5/home/epi/bhattacharya_lab/data/GTEx_v9/SQANTI3/*_filtered/*_classification.txt
#       SQANTI3 classifications for GTEx v9 tissues
#
# Output Files:
#   Figures:
#     - Figures/RIN.png
#         RIN score comparison between GDM and control groups for ONT and GUSTO cohorts
#     - Figures/fig1b.pdf
#         Bar plot comparing gene and isoform counts between GENCODEv45 and lr-assembly
#     - Figures/fig1c.pdf
#         Structural category distribution before and after filtering
#     - Figures/figs1a.png
#         Transcript length distribution by structural category
#     - Figures/figs1b.pdf
#         Exon count distribution by structural category
#     - Figures/figs1c.pdf
#         Forest plot of mean isoforms per gene comparing lr-assembly and GENCODEv45
#     - Figures/fig1d_revision_a.pdf
#         Transcriptional breadth (genes and isoforms) across GTEx tissues and placenta
#     - Figures/fig1e_revision_b.pdf
#         Saturation curves showing isoforms detected at varying expression thresholds
#     - Figures/fig1f_revision.pdf
#         Transcriptional complexity (isoforms per gene) across tissues
#
################################################################################

################################################################################
# Environment Setup
################################################################################
# .libPaths(c( "/home/stbresnahan/R/ubuntu/4.3.1" , .libPaths()))

# Load required libraries for data manipulation, genomics, and visualization
library(dplyr)         # Data manipulation
library(rtracklayer)   # Genomic data import/export
library(ggplot2)       # Data visualization
library(tidyr)         # Data tidying
library(cowplot)       # Plot arrangement
library(plyr)          # Data manipulation (legacy functions)
library(ggExtra)       # Additional ggplot2 functionality
library(matrixStats)   # Matrix statistics
library(biomaRt)       # Biomart database queries
library(stringr)       # String manipulation
library(scales)        # Scale functions for ggplot2
library(patchwork)     # For making multi-panel figures
library(data.table)    # For parallel operations on tables
library(parallel)      # For parallel processing

### RIN scores ###

# ONT data
agilent_data <- read.csv("Supplementary/ONT_RIN_scores.csv")
t_test_ont <- t.test(RIN ~ GDM, data = agilent_data)
label_ont <- paste0("Difference: ", round(diff(t_test_ont$estimate), 2), 
                    "\np = ", round(t_test_ont$p.value, 3))

g1 <- ggplot(agilent_data, aes(x = GDM, y = RIN)) + 
  geom_boxplot() +
  annotate("text", x = 1.5, y = max(agilent_data$RIN, na.rm = TRUE), 
           label = label_ont, vjust = 1) +
  theme_minimal() + labs(title = "Long-read RIN scores") +
  theme(axis.title.x=element_blank())
g1

# GUSTO data
GUSTO_data <- read.csv("Supplementary/GUSTO_RIN_scores.csv")
GUSTO_data <- GUSTO_data[!is.na(GUSTO_data$GDM.Status), ]
t_test_gusto <- t.test(RIN ~ GDM.Status, data = GUSTO_data)
label_gusto <- paste0("Difference: ", round(diff(t_test_gusto$estimate), 2), 
                      "\np = ", round(t_test_gusto$p.value, 3))

g2 <- ggplot(GUSTO_data, aes(x = GDM.Status, y = RIN)) + 
  geom_boxplot() +
  annotate("text", x = 1.5, y = max(GUSTO_data$RIN, na.rm = TRUE), 
           label = label_gusto, vjust = 1) +
  theme_minimal() + labs(title = "GUSTO short-read RIN scores") +
  theme(axis.title.x=element_blank())
g2

g3 <- g1 + g2 + plot_layout(ncol = 2)

ggsave("Figures/RIN.png",g3,width=6.57,height=3.85,units="in",dpi=300)


################################################################################
# Data Import and Processing
################################################################################

# Import SQANTI3 classification results
# SQANTI3 classifies transcripts based on their structural relationship to reference
classif <- read.table("Assembly/ESPRESSO_corrected_classification.txt", header = TRUE)

# Standardize column naming
names(classif)[1] <- "transcript_id"

# Define structural category hierarchy and factor levels
# Categories represent increasing deviation from reference annotation
classif$structural_category <- factor(as.character(classif$structural_category),
                                      levels = c("full-splice_match",
                                                 "incomplete-splice_match", 
                                                 "novel_in_catalog",
                                                 "novel_not_in_catalog",
                                                 "genic",
                                                 "antisense", 
                                                 "fusion",
                                                 "intergenic",
                                                 "genic_intron"))

# Rename categories for cleaner visualization
# FSM: Full Splice Match - exact match to reference
# ISM: Incomplete Splice Match - subset of reference exons
# NIC: Novel In Catalog - novel combination of known splice sites
# NNC: Novel Not in Catalog - contains novel splice sites
classif$structural_category <- revalue(classif$structural_category, 
                                       c("full-splice_match" = "FSM",
                                         "incomplete-splice_match" = "ISM",
                                         "novel_in_catalog" = "NIC", 
                                         "novel_not_in_catalog" = "NNC",
                                         "genic" = "Genic Genomic",
                                         "antisense" = "Antisense",
                                         "fusion" = "Fusion",
                                         "intergenic" = "Intergenic",
                                         "genic_intron" = "Genic Intron"))

# Identify transcripts with open reading frames (ORFs)
classif$ORF <- FALSE
classif[!is.na(classif$ORF_seq), "ORF"] <- TRUE

# Define color palette for structural categories (ColorBrewer Set3)
myPalette <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#ffd92f",
               "#e5c494", "#c87774", "#d3adaf", "#b3b3b3")

################################################################################
# Transcript-to-Gene Mapping for Long-Read Assembly
################################################################################

# Import GTF annotation for the corrected long-read assembly
tx2g.assembly <- data.frame(import("Assembly/ESPRESSO_corrected_corrected.gtf"))

# Clean and standardize the transcript-to-gene mapping
tx2g.assembly <- tx2g.assembly[!duplicated(tx2g.assembly), ]
tx2g.assembly$score <- NULL
tx2g.assembly$phase <- NULL

# Handle novel transcripts without gene assignments
novelTx.assembly <- unique(tx2g.assembly[tx2g.assembly$gene_id == "NA", "transcript_id"])

# Use SQANTI3 associated_gene assignments for novel transcripts when available
novelTx.ascg <- classif[classif$transcript_id %in% novelTx.assembly, 
                        c("transcript_id", "associated_gene")]
# Keep only those with valid ENSEMBL gene IDs
novelTx.ascg <- novelTx.ascg[grep("ENSG", novelTx.ascg$associated_gene), ]

# Merge associated gene information
tx2g.assembly <- left_join(tx2g.assembly, novelTx.ascg, by = "transcript_id")
tx2g.assembly[is.na(tx2g.assembly$associated_gene), "associated_gene"] <- "NA"
tx2g.assembly[tx2g.assembly$gene_id == "NA", "gene_id"] <- 
  tx2g.assembly[tx2g.assembly$gene_id == "NA", "associated_gene"]
tx2g.assembly$associated_gene <- NULL

# Create gene IDs for remaining unannotated novel transcripts
# Use genomic coordinates (chr:start:end) as gene identifier
novelTx.assembly.unannot <- setdiff(novelTx.assembly, novelTx.ascg$transcript_id)
tx2g.assembly[tx2g.assembly$transcript_id %in% novelTx.assembly.unannot, "gene_id"] <- 
  sapply(tx2g.assembly[tx2g.assembly$transcript_id %in% novelTx.assembly.unannot, "transcript_id"],
         function(x) paste(strsplit(x, ":")[[1]][1:3], collapse = ":"))

# Finalize transcript-to-gene mapping
tx2g.assembly <- tx2g.assembly[, c("transcript_id", "gene_id")]
tx2g.assembly <- tx2g.assembly[!duplicated(tx2g.assembly), ]

################################################################################
# Transcript-to-Gene Mapping for GENCODE v45 Reference
################################################################################

# Import GENCODE v45 reference annotation
tx2g.gencode <- data.frame(import("/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/GenomicReferences/txome/gencode_v45/gencode.v45.annotation.gtf"))

# Filter for standard chromosomes only (exclude scaffolds and patches)
tx2g.gencode <- tx2g.gencode[tx2g.gencode$seqnames %in% 
                               c(paste0("chr", seq.int(1, 22, 1)), "chrX", "chrY", "chrM"), ]

# Extract transcript-to-gene mapping
tx2g.gencode <- tx2g.gencode[, c("transcript_id", "gene_id")]
tx2g.gencode <- tx2g.gencode[!is.na(tx2g.gencode$transcript_id), ]
tx2g.gencode <- tx2g.gencode[!duplicated(tx2g.gencode), ]

################################################################################
# Figure 1B: Gene and Isoform Counts - Assembly vs Reference Comparison
################################################################################

# Calculate counts for visualization
fig1b.dat <- data.frame(
  Annotation = c("GENCODEv45", "GENCODEv45", "lr-assembly", "lr-assembly"),
  Set = factor(c("Genes", "Isoforms", "Genes", "Isoforms"), levels = c("Isoforms", "Genes")),
  Count = c(
    length(unique(tx2g.gencode$gene_id)),      # GENCODE genes
    length(unique(tx2g.gencode$transcript_id)), # GENCODE transcripts
    length(unique(tx2g.assembly$gene_id)),     # Assembly genes  
    length(unique(tx2g.assembly$transcript_id)) # Assembly transcripts
  ),
  # Calculate percentage reduction from reference
  Labels = c("", "", "-63.5%", "-73.1%")
)
# Create bar plot comparing gene and isoform counts
fig1b <- ggplot(fig1b.dat, aes(y = Count, x = Annotation, fill = Annotation, group = Set)) +
  geom_bar(stat = "identity", position = "dodge", aes(alpha = Set)) + 
  scale_fill_manual(values = c("GENCODEv45" = "royalblue", "lr-assembly" = "seagreen")) +
  scale_alpha_manual(values = c("Isoforms" = 0.5, "Genes" = 1),
                     guide = guide_legend(override.aes = list(fill = "grey30"))) +
  theme_minimal() + 
  scale_y_continuous(labels = scales::label_number(scale = 1e-3, suffix = "k")) +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = Labels), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, 
            color = "black",
            show.legend = FALSE) +
  guides(fill = "none") + ggtitle("Annotated features")
fig1b

ggsave("Figures/fig1b.pdf", fig1b, width = 3.56, height = 4.32, units = "in", dpi = 300)

################################################################################
# Figure 1C: Structural Category Distribution
################################################################################

# Count transcripts by structural category
fig1c.dat.a <- data.frame(table(classif$structural_category))
names(fig1c.dat.a)[1] <- "structural_category"

# Count transcripts by structural category after filtering
classif.filtered <- read.table("Assembly/ESPRESSO_corrected_SC_filtered_classification.txt", header = TRUE)
names(classif.filtered)[1] <- "transcript_id"
classif.filtered$structural_category <- factor(as.character(classif.filtered$structural_category),
                                      levels = c("full-splice_match",
                                                 "incomplete-splice_match", 
                                                 "novel_in_catalog",
                                                 "novel_not_in_catalog",
                                                 "genic",
                                                 "antisense", 
                                                 "fusion",
                                                 "intergenic",
                                                 "genic_intron"))
classif.filtered$structural_category <- revalue(classif.filtered$structural_category, 
                                       c("full-splice_match" = "FSM",
                                         "incomplete-splice_match" = "ISM",
                                         "novel_in_catalog" = "NIC", 
                                         "novel_not_in_catalog" = "NNC",
                                         "genic" = "Genic Genomic",
                                         "antisense" = "Antisense",
                                         "fusion" = "Fusion",
                                         "intergenic" = "Intergenic",
                                         "genic_intron" = "Genic Intron"))

fig1c.dat.b <- data.frame(table(classif.filtered$structural_category))
names(fig1c.dat.b)[1] <- "structural_category"

# Create bar plot of structural categories
fig1c.dat <- rbind(
  transform(fig1c.dat.a, Group = "Before"),
  transform(fig1c.dat.b, Group = "After")
)
fig1c.dat$Main <- ifelse(fig1c.dat$structural_category %in% c("FSM", "ISM", "NIC", "NNC"), 
                         "Main", "Other")
fig1c.dat$Group <- factor(fig1c.dat$Group, levels = c("Before", "After"))

fig1c <- ggplot(fig1c.dat, aes(x = structural_category, y = Freq, 
                               fill = structural_category, alpha = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = myPalette) +
  scale_alpha_manual(values = c("Before" = 0.5, "After" = 1),
                     labels = c("All Isoforms", "High-quality isoforms"),
                     name = NULL) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-3, suffix = "k")) +
  theme_minimal() +
  theme(legend.position = "top",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        text = element_text(color = "black"),
        axis.text = element_text(color = "black")) +
  ylab("Count") +
  ggtitle("Isoform categories") +
  facet_wrap(~ Main, scales = "free") +
  guides(fill = "none")
fig1c
ggsave("Figures/fig1c.pdf", fig1c, width = 4.74, height = 4.32, units = "in", dpi = 300)

################################################################################
# Figure S1A: Transcript Length Distribution by Structural Category
################################################################################

# Box plot of transcript lengths by structural category
figs1a <- ggplot(classif, aes(y = length, x = structural_category, fill = structural_category)) +
  geom_boxplot() + 
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_fill_manual(values = myPalette) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "none") +
  labs(title = "Isoform length", y = "Length (bp)")

figs1a
ggsave("Figures/figs1a.png", figs1a, width = 3, height = 3, units = "in", dpi = 300)

################################################################################
# Figure S1B: Exon Count Distribution by Structural Category  
################################################################################

# Stacked histogram of exon counts by structural category
figs1b <- ggplot(classif, aes(x = structural_category, y = exons, fill = structural_category)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # suppress outliers if needed
  scale_y_continuous(limits = c(0, 25), labels = scales::comma) +
  scale_fill_manual(values = myPalette) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  labs(title = "Exons per isoform by structural category",
       x = "Structural category", 
       y = "Number of exons")

figs1b

ggsave("Figures/figs1b.pdf", figs1b, width = 4, height = 4, units = "in", dpi = 300)

################################################################################
# Figure S1C: Isoform Complexity Comparison - Assembly vs Reference
################################################################################

# Calculate isoforms per gene for assembly
iso.tab.assembly <- data.frame(table(tx2g.assembly$gene_id))
names(iso.tab.assembly) <- c("gene_id", "Freq.assembly")
iso.tab.assembly$gene_id <- as.character(iso.tab.assembly$gene_id)

# Calculate isoforms per gene for GENCODE reference
iso.tab.gencode <- data.frame(table(tx2g.gencode$gene_id))
names(iso.tab.gencode) <- c("gene_id", "Freq.gencode")
iso.tab.gencode$gene_id <- as.character(iso.tab.gencode$gene_id)

# Create comprehensive comparison dataset
genes.all <- union(tx2g.assembly$gene_id, tx2g.gencode$gene_id)
compgene <- data.frame(gene_id = genes.all)
compgene <- left_join(compgene, iso.tab.assembly, by = "gene_id")
compgene <- left_join(compgene, iso.tab.gencode, by = "gene_id")
compgene[is.na(compgene)] <- 0

# Classify genes by isoform complexity change
compgene$set <- "no change"
compgene[(compgene$Freq.gencode - compgene$Freq.assembly) > 0, "set"] <- "reduced"
compgene[(compgene$Freq.gencode - compgene$Freq.assembly) < 0, "set"] <- "expanded"
compgene$set <- factor(compgene$set, levels = c("expanded", "reduced", "no change"))
compgene <- compgene[order(compgene$set), ]

# Pivot to long format for boxplot
compgene_long <- compgene %>%
  dplyr::select(Freq.assembly, Freq.gencode) %>%
  pivot_longer(cols = everything(),
               names_to = "dataset",
               values_to = "isoforms") %>%
  mutate(dataset = recode(dataset,
                          "Freq.assembly" = "lr-assembly",
                          "Freq.gencode" = "GENCODEv45"))

# Summarize mean, SD, and max per dataset
compgene_summary <- compgene_long %>%
  mutate(log_isoforms = log2(isoforms + 1)) %>%
  group_by(dataset) %>%
  summarise(
    mean_val = mean(log_isoforms),
    sd_val = sd(log_isoforms),
    max_val = max(isoforms)
  )

compgene_summary <- compgene_long %>%
  group_by(dataset) %>%
  summarise(
    mean_val = mean(isoforms),
    sd_val = sd(isoforms),
    max_val = max(isoforms)
  )

myColors <- c("lr-assembly" = "seagreen", "GENCODEv45" = "royalblue")

# Forest-style plot
figs1c <- ggplot(compgene_summary, aes(y = dataset, x = mean_val, color = dataset)) +
  geom_point(size = 4) +
  geom_errorbarh(aes(xmin = mean_val - sd_val, xmax = mean_val + sd_val), height = 0.2) +
  geom_text(aes(x = mean_val + sd_val + 0.2, label = paste0("(", max_val, ")")), 
            hjust = 0, size = 3,color="black") +
  scale_color_manual(values = myColors) +
  theme_minimal() +
  xlim(-4,14) +
  labs(x = "Mean count of isoforms per gene ± SD",
       y = "",
       title = "Isoforms per gene by dataset")

figs1c

# Save as PDF
ggsave("Figures/figs1c.pdf", 
       plot = figs1c,
       width = 4.8, 
       height = 4.8, 
       units = "in",
       dpi = 300,
       device = "pdf")

############################################################################################################
# Figure 1D: Transcriptional breadth across tissues
############################################################################################################

# Get all directories in the GTEx folder
gtex_path <- "/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v9/SQANTI3"
tissues <- list.dirs(gtex_path, full.names = FALSE, recursive = FALSE)
tissues <- tissues[grep("filtered",tissues)]

# Initialize result lists
isoforms_per_tissue <- list()
genes_per_tissue <- list()
gene_isoform_counts <- list()

# Loop through each tissue directory
for (tissue in tissues) {
  # Construct the classification file path
  classif_file <- file.path(gtex_path, tissue, 
                            paste0(tissue, "_classification.txt"))
  
  # Check if file exists
  if (file.exists(classif_file)) {
    # Read classification table
    classif <- read.table(classif_file, header = TRUE, sep = "\t")
    
    # Count unique isoforms and genes
    isoforms_per_tissue[[tissue]] <- length(unique(classif$isoform))
    genes_per_tissue[[tissue]] <- length(unique(classif$associated_gene))
    
    # Count isoforms per gene
    gene_counts <- aggregate(isoform ~ associated_gene, 
                             data = classif, 
                             FUN = function(x) length(unique(x)))
    colnames(gene_counts) <- c("gene", "isoforms")
    gene_counts$tissue <- tissue
    gene_isoform_counts[[tissue]] <- gene_counts[, c("tissue", "gene", "isoforms")]
  }
}

# Combine gene-isoform counts into a single data.frame
gene_isoform_df <- do.call(rbind, gene_isoform_counts)
rownames(gene_isoform_df) <- NULL

# Create final result list
result <- list(
  isoforms_per_tissue = isoforms_per_tissue,
  genes_per_tissue = genes_per_tissue,
  gene_isoform_counts = gene_isoform_df
)

counts_df <- data.frame(
  tissue = names(result$isoforms_per_tissue),
  genes = unlist(result$genes_per_tissue),
  isoforms = unlist(result$isoforms_per_tissue),
  row.names = NULL
)

dat.fig1d <- reshape(counts_df, 
                     direction = "long",
                     varying = c("genes", "isoforms"),
                     v.names = "counts",
                     timevar = "level",
                     times = c("genes", "isoforms"),
                     idvar = "tissue")
rownames(dat.fig1d) <- NULL
dat.fig1d <- data.frame(dat.fig1d)
dat.fig1d[dat.fig1d$level == "genes", "level"] <- "Genes"
dat.fig1d[dat.fig1d$level == "isoforms", "level"] <- "Isoforms"

dat.fig1d$tissue <- dat.fig1d$tissue %>%
  str_remove("_filtered") %>%
  str_replace_all("_", " ") %>%
  str_replace("(Brain|Adipose|Breast|Cells) ", "\\1 - ") %>%
  str_replace(" (BA\\d+)", " (\\1)") %>%
  str_to_title()

dat.fig1d <- rbind(dat.fig1d,data.frame(tissue=c("Placenta","Placenta"),
                                        level=c("Genes","Isoforms"),
                                        counts=c(12302,37661)))

dat.fig1d$tissue <- c(paste0(dat.fig1d$tissue[1]," (n = 1.8M)"),
                      paste0(dat.fig1d$tissue[2]," (n = 2.8M)"),
                      paste0(dat.fig1d$tissue[3]," (n = 3.5M)"),
                      paste0(dat.fig1d$tissue[4]," (n = 48.5M)"),
                      paste0(dat.fig1d$tissue[5]," (n = 22.9M)"),
                      paste0(dat.fig1d$tissue[6]," (n = 20.2M)"),
                      paste0(dat.fig1d$tissue[7]," (n = 2.3M)"),
                      paste0(dat.fig1d$tissue[8]," (n = 83.9M)"),
                      paste0(dat.fig1d$tissue[9]," (n = 41.3M)"),
                      paste0(dat.fig1d$tissue[10]," (n = 20.2M)"),
                      paste0(dat.fig1d$tissue[11]," (n = 7.6M)"),
                      paste0(dat.fig1d$tissue[12]," (n = 37.3M)"),
                      paste0(dat.fig1d$tissue[13]," (n = 41.5M)"),
                      paste0(dat.fig1d$tissue[14]," (n = 41.0M)"),
                      paste0(dat.fig1d$tissue[15]," (n = 6.0M)"),
                      paste0(dat.fig1d$tissue[1]," (n = 1.8M)"),
                      paste0(dat.fig1d$tissue[2]," (n = 2.8M)"),
                      paste0(dat.fig1d$tissue[3]," (n = 3.5M)"),
                      paste0(dat.fig1d$tissue[4]," (n = 48.5M)"),
                      paste0(dat.fig1d$tissue[5]," (n = 22.9M)"),
                      paste0(dat.fig1d$tissue[6]," (n = 20.2M)"),
                      paste0(dat.fig1d$tissue[7]," (n = 2.3M)"),
                      paste0(dat.fig1d$tissue[8]," (n = 83.9M)"),
                      paste0(dat.fig1d$tissue[9]," (n = 41.3M)"),
                      paste0(dat.fig1d$tissue[10]," (n = 20.2M)"),
                      paste0(dat.fig1d$tissue[11]," (n = 7.6M)"),
                      paste0(dat.fig1d$tissue[12]," (n = 37.3M)"),
                      paste0(dat.fig1d$tissue[13]," (n = 41.5M)"),
                      paste0(dat.fig1d$tissue[14]," (n = 41.0M)"),
                      paste0(dat.fig1d$tissue[15]," (n = 6.0M)"),
                      paste0(dat.fig1d$tissue[31]," (n = 118.6M)"),
                      paste0(dat.fig1d$tissue[32]," (n = 118.6M)"))
dat.fig1d$tissue <- factor(dat.fig1d$tissue, levels = sort(unique(dat.fig1d$tissue)))

# Plot
n_tissues <- length(unique(dat.fig1d$tissue))

dat.fig1d <- dat.fig1d %>%
  mutate(
    n_label = str_extract(tissue, "(?<=\\(n = )([0-9\\.]+)([MK]?)"),
    n_numeric = as.numeric(str_extract(n_label, "[0-9\\.]+")),
    n_multiplier = case_when(
      str_detect(n_label, "M") ~ 1e6,
      str_detect(n_label, "K") ~ 1e3,
      TRUE ~ 1
    ),
    n_value = n_numeric * n_multiplier
  ) %>%
  arrange(desc(n_value))  # sort highest to lowest

tissue_levels <- dat.fig1d %>%
  arrange(desc(n_value)) %>%
  pull(tissue) %>%
  unique()

dat.fig1d <- dat.fig1d %>%
  mutate(
    tissue = factor(tissue, levels = rev(tissue_levels))
  )

tissue_colors <- setNames(
  scales::hue_pal()(nrow(dat.fig1d)),
  dat.fig1d$tissue
)

tissue_colors["Placenta (n = 118.6M)"] <- "black"

fig1d <- ggplot(dat.fig1d, aes(x = tissue, y = counts, fill = tissue, alpha = level)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(
    title = "Transcriptional breadth across tissues",
    subtitle = "n = assembled reads",
    y = "Count"
  ) +
  scale_y_continuous(limits = c(0, 50000),
                     labels = label_number(scale = 1e-3, suffix = "k")) +
  scale_fill_manual(values = tissue_colors, guide = "none") +
  scale_alpha_manual(values = c(0.5, 1), 
                     breaks = c("Isoforms", "Genes"),
                     guide = guide_legend(override.aes = list(fill = "grey30"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        legend.position = c(0.8, 0.1),
        legend.title = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(hjust = 0.75),
        plot.subtitle = element_text(hjust = 0.75)) + 
  coord_flip()
fig1d

ggsave("Figures/fig1d_revision_a.pdf", fig1d, width = 8.42, height = 3.74, units = "in", dpi = 300)


############################################################################################################
# Figure 1E: Saturation curves for read depth to transcripts detected
############################################################################################################

# Function to process each tissue
process_tissue <- function(tissue, gtex_path) {
  classif_file <- file.path(gtex_path, tissue, 
                            paste0(tissue, "_classification.txt"))
  
  if (file.exists(classif_file)) {
    classif <- read.table(classif_file, header = TRUE, sep = "\t")
    
    # Count isoforms at each FL value
    isoform_counts <- table(classif$FL)
    fl_values <- as.numeric(names(isoform_counts))
    
    # Create bins of 10
    max_fl <- max(fl_values)
    fl_thresholds <- seq(10, ceiling(max_fl / 10) * 10, by = 10)
    
    # Cumulative count of isoforms at each bin
    results <- data.frame(
      tissue = tissue,
      FL_threshold = fl_thresholds,
      Isoforms = sapply(fl_thresholds, function(fl) sum(isoform_counts[fl_values <= fl]))
    )
    
    return(results)
  }
  return(NULL)
}

# Process all tissues in parallel
n_cores <- detectCores() - 1
saturation_data <- mclapply(tissues, process_tissue, gtex_path = gtex_path, mc.cores = n_cores)
saturation_data <- do.call(rbind, saturation_data)

# Add Placenta from classif.filtered
isoform_counts_placenta <- table(classif.filtered$FL)
fl_values_placenta <- as.numeric(names(isoform_counts_placenta))
max_fl_placenta <- max(fl_values_placenta)
fl_thresholds_placenta <- seq(10, ceiling(max_fl_placenta / 10) * 10, by = 10)

placenta_data <- data.frame(
  tissue = "Placenta",
  FL_threshold = fl_thresholds_placenta,
  Isoforms = sapply(fl_thresholds_placenta, function(fl) sum(isoform_counts_placenta[fl_values_placenta <= fl]))
)

# Combine into single data frame
sat_df <- rbind(saturation_data, placenta_data)

# Clean tissue names
sat_df$tissue <- sat_df$tissue %>%
  str_remove("_filtered") %>%
  str_replace_all("_", " ") %>%
  str_replace("(Brain|Adipose|Breast|Cells) ", "\\1 - ") %>%
  str_replace(" (BA\\d+)", " (\\1)") %>%
  str_to_title()

# Remove (n = …) from names and keep first occurrence
tissue_names <- gsub(" \\(n = .*\\)", "", names(tissue_colors))
tissue_colors_sat <- tissue_colors[!duplicated(tissue_names)]  # keep first
names(tissue_colors_sat) <- tissue_names[!duplicated(tissue_names)]

# Make sat_df$tissue a factor with these unique levels
sat_df <- sat_df %>%
  mutate(tissue = factor(tissue, levels = names(tissue_colors_sat)))

# Linewidth mapping
linewidths <- rep(0.5, length(names(tissue_colors_sat)))
names(linewidths) <- names(tissue_colors_sat)
linewidths["Placenta"] <- 1.5

# Plot
fig1e <- ggplot(sat_df, aes(x = FL_threshold, y = Isoforms,
                              color = tissue, linewidth = tissue)) +
  geom_line() +
  labs(
    x = "Expression",
    y = "Count"
  ) +
  scale_x_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  scale_color_manual(values = tissue_colors_sat) +
  scale_linewidth_manual(values = linewidths) +
  theme_minimal() +
  ggtitle("Isoforms detected at expression threshold") +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    legend.title = element_blank(),
    legend.position = "none"
  )

fig1e

ggsave("Figures/fig1e_revision_b.pdf", fig1d_b, width = 4.5, height = 3.82, units = "in", dpi = 300)

################################################################################
# Figure 1F: Transcriptional Complexity Across Tissues
################################################################################

dat.fig1e <- result[["gene_isoform_counts"]]

# Clean tissue names
dat.fig1e$tissue <- dat.fig1e$tissue %>%
  str_remove("_filtered") %>%
  str_replace_all("_", " ") %>%
  str_replace("(Brain|Adipose|Breast|Cells) ", "\\1 - ") %>%
  str_replace(" (BA\\d+)", " (\\1)") %>%
  str_to_title()

# Add placenta
dat.fig1e.placenta <- data.frame(table(classif.filtered$associated_gene))
dat.fig1e.placenta$tissue <- "Placenta"
names(dat.fig1e.placenta)[1:2] <- c("gene","isoforms") 
dat.fig1e.placenta <- dat.fig1e.placenta[,c(3,1:2)]
dat.fig1e <- rbind(dat.fig1e, dat.fig1e.placenta)

# ----------------------------
# Use tissue order from saturation curves (sat_df)
# ----------------------------
tissue_levels <- c("Placenta",levels(sat_df$tissue)[-grep("Placenta",levels(sat_df$tissue))])
dat.fig1e$tissue <- factor(dat.fig1e$tissue, levels = rev(tissue_levels))  # reverse for coord_flip

# ----------------------------
# Use the same color palette as Figure 1D
# ----------------------------
palette_colors <- tissue_colors_sat[tissue_levels]
palette_colors <- palette_colors[names(palette_colors) %in% tissue_levels]  # keep only relevant tissues

# ----------------------------
# Compute range labels for max isoform per tissue
# ----------------------------
range_labels <- dat.fig1e %>%
  dplyr::group_by(tissue) %>%
  dplyr::summarise(max_val = max(as.numeric(isoforms), na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    label = paste0("(", max_val, ")"),
    # place label slightly above the SD line
    y_pos = 9
  )

# ----------------------------
# Summarize mean and SD per tissue, clamp lower bound at 0
# ----------------------------
dat_summary <- dat.fig1e %>%
  dplyr::group_by(tissue) %>%
  dplyr::summarise(
    mean_isoforms = mean(as.numeric(isoforms), na.rm = TRUE),
    sd_isoforms = sd(as.numeric(isoforms), na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    ymin = pmax(mean_isoforms - sd_isoforms, 0),  # clamp lower bound at 0
    ymax = mean_isoforms + sd_isoforms
  )

# ----------------------------
# Plot mean ± SD with max labels
# ----------------------------

fig1f <- ggplot(dat_summary, aes(x = tissue, y = mean_isoforms, color = tissue)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2) +
  geom_text(data = range_labels, aes(x = tissue, y = y_pos, label = label),
            hjust = 0, size = 3, inherit.aes = FALSE) +
  scale_color_manual(values = palette_colors) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_text(color = "black"),
    plot.title.position = "plot",   # left-align title
    plot.title = element_text(hjust = 0),
    plot.subtitle = element_text(hjust = 0)
  ) +
  coord_flip() +
  scale_y_continuous(limits = c(0,10), breaks = seq(0, 10, 2)) +
  labs(
    y = "Isoforms per gene",
    title = "Transcriptional complexity across tissues",
    subtitle = "Points = mean; lines = ±1 SD; (max) shown above"
  )

fig1f


ggsave("Figures/fig1f_revision.pdf", fig1e, width = 4.58, height = 3.74, units = "in", dpi = 300)

################################################################################
# End of Analysis
################################################################################