################################################################################
# Placental Long-Read RNA-seq Transcriptome Assembly Analysis
# Author: Sean T. Bresnahan
# Description: Analysis of placental transcriptome assembly against GENCODEv45
#              reference using SQANTI3 classification and comparison with GTEx v9
#              
# Input Files:
#   - SQANTI3 classification results (ESPRESSO_corrected_classification.txt)
#   - Raw assembled GTF (ESPRESSO_corrected_corrected.gtf)
#   - GENCODE v45 reference annotation (gencode.v45.annotation.gtf)
#   - GTEx v9 expression data (quantification_gencode.counts.txt)
#   - GTEx v9 sample metadata (GTEx_Sample_Attributes.GRU.txt)
#
# Output Figures:
#   - Main Manuscript: 1B-1E
#     * 1B: Gene and isoform counts (assembly vs reference)
#     * 1C: Structural category distribution
#     * 1D: Multi-tissue gene and isoform detection comparison
#     * 1E: Transcriptional complexity across tissues
#   - Supplementary: S1A-S1C
#     * S1A: Transcript length distribution by structural category
#     * S1B: Exon count distribution by structural category
#     * S1C: Isoform complexity comparison (assembly vs reference)
################################################################################

################################################################################
# Environment Setup
################################################################################

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

################################################################################
# Data Import and Processing
################################################################################

# Import SQANTI3 classification results
# SQANTI3 classifies transcripts based on their structural relationship to reference
classif <- read.table("ESPRESSO_corrected_classification.txt", header = TRUE)

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
tx2g.assembly <- data.frame(import("ESPRESSO_corrected_corrected.gtf"))

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
tx2g.gencode <- data.frame(import("gencode.v45.annotation.gtf"))

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
  Set = c("Genes", "Isoforms", "Genes", "Isoforms"),
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
fig1b <- ggplot(fig1b.dat, aes(y = Count, x = Annotation, fill = Set)) +
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_manual(values = c("black", "grey70")) +
  theme_minimal() + 
  scale_y_continuous(labels = scales::comma) +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        axis.title.x = element_blank()) +
  geom_text(aes(label = Labels), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, 
            color = "black")

# Display and save figure
fig1b
ggsave("fig1b.png", fig1b, width = 3.5, height = 4, units = "in", dpi = 300)

################################################################################
# Figure 1C: Structural Category Distribution
################################################################################

# Count transcripts by structural category
fig1c.dat <- data.frame(table(classif$structural_category))
names(fig1c.dat)[1] <- "structural_category"

# Create bar plot of structural categories
fig1c <- ggplot(fig1c.dat, aes(x = structural_category, y = Freq, fill = structural_category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = myPalette) +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  ylab("Isoforms") +
  ggtitle("Structural categories")

fig1c
ggsave("fig1c.png", fig1c, width = 3, height = 3, units = "in", dpi = 300)

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
ggsave("figs1a.png", figs1a, width = 3, height = 3, units = "in", dpi = 300)

################################################################################
# Figure S1B: Exon Count Distribution by Structural Category  
################################################################################

# Stacked histogram of exon counts by structural category
figs1b <- ggplot(classif, aes(x = exons, fill = structural_category)) +
  geom_bar(stat = "count", position = "stack") +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values = myPalette) +
  theme_minimal() + 
  theme(legend.position = "none") +
  labs(title = "Exons per isoform by structural category",
       x = "Exons", 
       y = "Isoforms") +
  xlim(0, 25)  # Focus on reasonable exon count range

figs1b
ggsave("figs1b.png", figs1b, width = 4, height = 4, units = "in", dpi = 300)

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

# Highlight specific genes of interest (Growth Hormone genes)
compgene$col <- ifelse(compgene$gene_id %in% c("ENSG00000136488.15", "ENSG00000136487.18"), 
                       "highlight", "normal")
highlighted <- compgene[compgene$col == "highlight", ]

# Create scatter plot with marginal histograms
g6 <- ggplot(compgene, aes(x = log(Freq.gencode + 1), y = log(Freq.assembly + 1))) +
  geom_point(size = 0.6, aes(color = col)) +
  geom_label(data = highlighted,
             aes(label = c("GH2", "CSH1")),
             color = "black",
             fill = "white", 
             label.size = 0.3,
             size = 3,
             label.r = unit(0.1, "lines"),
             nudge_y = 0.25) +
  scale_color_manual(values = c("highlight" = "red", "normal" = "grey60")) +
  xlim(0, 6) + ylim(0, 6) +
  geom_abline(slope = 1, intercept = 0, color = 'black',
              linetype = "dashed", linewidth = 0.5) +
  theme_minimal() + 
  theme(legend.position = "none") +
  xlab("GENCODEv45 (x10³)") +
  ylab("lr-assembly (x10³)") +
  ggtitle("Isoforms per gene")

# Add marginal histograms
ggMarginal(g6, type = "histogram", fill = "grey40", color = "white")

################################################################################
# GTEx v9 Multi-Tissue Comparison Analysis
################################################################################

# Import GTEx v9 expression data and sample metadata
GTEx_counts <- read.table("Supplementary/quantification_gencode.counts.txt", header = TRUE)
names(GTEx_counts) <- gsub("[.]", "-", names(GTEx_counts))

GTEx_meta <- read.table("Supplementary/GTEx_Sample_Attributes.GRU.txt", 
                        header = TRUE, sep = "\t", fill = NA)
GTEx_meta <- GTEx_meta[, c("SAMPID", "SMTSD")]
names(GTEx_meta) <- c("SampleID", "Tissue")

# Aggregate expression data by tissue type
transcript_ids <- GTEx_counts[[1]]
expr_matrix <- as.matrix(GTEx_counts[, -1])
colnames(expr_matrix) <- colnames(GTEx_counts)[-1]

# Create sample-to-tissue mapping
sample_to_tissue <- setNames(GTEx_meta$Tissue, GTEx_meta$SampleID)
common_samples <- intersect(colnames(expr_matrix), names(sample_to_tissue))
expr_matrix <- expr_matrix[, common_samples]
tissue_labels <- sample_to_tissue[common_samples]

# Calculate mean expression per tissue
tissues <- unique(tissue_labels)
mean_by_tissue <- sapply(tissues, function(tiss) {
  sample_cols <- which(tissue_labels == tiss)
  rowMeans(expr_matrix[, sample_cols, drop = FALSE], na.rm = TRUE)
})

GTEx_mean_by_tissue <- data.frame(transcript_id = transcript_ids, 
                                  mean_by_tissue, 
                                  check.names = FALSE)

# Handle missing values and apply expression threshold
GTEx_mean_by_tissue <- GTEx_mean_by_tissue %>%
  mutate(across(-transcript_id, ~replace_na(., 0)))
GTEx_mean_by_tissue[is.na(GTEx_mean_by_tissue)] <- 0

# Convert to long format and filter by expression threshold (>3 counts)
long_df <- as.data.frame(GTEx_mean_by_tissue) %>%
  pivot_longer(
    cols = -transcript_id,
    names_to = "Tissue",
    values_to = "Expression",
    values_transform = list(Expression = as.numeric, Tissue = as.character)
  )
long_df <- long_df[long_df$Expression > 3, ]

################################################################################
# Transcript-to-Gene Mapping for GTEx Analysis
################################################################################

# Query Ensembl for transcript-to-gene mapping
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
transcript_ids <- as.character(data.frame(unique(long_df$transcript_id)))
transcript_ids <- unlist(strsplit(gsub('[\"\\n]', '', transcript_ids), ",\\s*"))

# Remove version numbers for Ensembl query
transcript_ids_clean <- gsub("\\.[0-9]+$", "", transcript_ids)

# Get transcript-to-gene mapping from Ensembl
mapping <- getBM(
  attributes = c("ensembl_transcript_id", "ensembl_gene_id"),
  filters = "ensembl_transcript_id", 
  values = transcript_ids_clean,
  mart = ensembl
)
names(mapping) <- c("transcript_id", "gene_id")

# Process GTEx data for gene-level analysis
long_df <- long_df %>%
  mutate(transcript_id_base = gsub("\\.[0-9]+$", "", transcript_id))
long_df <- long_df[, c(4, 2, 3)]
names(long_df)[1] <- "transcript_id"

################################################################################
# Figure 1D: Multi-Tissue Gene and Isoform Detection Comparison
################################################################################

# Handle transcript ID version matching
GTEx_mean_by_tissue$transcript_id <- str_remove(GTEx_mean_by_tissue$transcript_id, "\\.\\d+$")

# Identify transcripts with multiple gene mappings (exclude for clean analysis)
mapping.tab <- data.frame(table(mapping$transcript_id))
names(mapping.tab) <- c("transcript_id", "Freq")

# Join expression data with gene mapping
expression_with_gene <- GTEx_mean_by_tissue %>%
  left_join(mapping[!mapping$transcript_id %in% mapping.tab[mapping.tab$Freq > 1, "transcript_id"], ], 
            by = "transcript_id")

# Convert to long format and filter by expression
expression_long <- expression_with_gene %>%
  pivot_longer(
    cols = -c(transcript_id, gene_id),
    names_to = "tissue",
    values_to = "expression"
  )

expression_filtered <- expression_long %>%
  filter(expression > 3)
expression_filtered <- expression_filtered[!is.na(expression_filtered$gene_id), ]

# Calculate tissue-specific statistics
tissues <- unique(expression_filtered$tissue)
txbyg <- data.frame(matrix(ncol = 5, nrow = 0))
names(txbyg) <- c("tissue", "transcripts", "genes", "mean.txpg", "sd.txpg")

for(t in tissues) {
  df <- data.frame(table(expression_filtered[expression_filtered$tissue == t, "gene_id"]))
  txbyg <- rbind(txbyg, data.frame(
    tissue = t,
    transcripts = sum(df$Freq),
    genes = length(df$Var1),
    mean.txpg = mean(df$Freq),
    sd.txpg = sd(df$Freq),
    range.txpg.low = range(df$Freq)[1],
    range.txpg.high = range(df$Freq)[2]
  ))
}

# Add placenta data from current study
txbyg <- rbind(txbyg, data.frame(
  tissue = "Placenta",
  transcripts = 37661,
  genes = 12302, 
  mean.txpg = 3.06,
  sd.txpg = 3.43,
  range.txpg.low = 1,
  range.txpg.high = 108
))

# Prepare data for visualization
txbyg_long <- txbyg %>%
  dplyr::select(tissue, transcripts, genes) %>%
  tidyr::pivot_longer(cols = c(transcripts, genes), names_to = "type", values_to = "count")

txbyg_long <- data.frame(txbyg_long)
txbyg_long[txbyg_long$type == "genes", "type"] <- "Genes"
txbyg_long[txbyg_long$type == "transcripts", "type"] <- "Isoforms"

# Create tissue comparison plot
fig1d <- ggplot(txbyg_long, aes(x = tissue, y = count, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(
    title = "lr-defined genes and isoforms",
    y = "Count",
    fill = "Type"
  ) +
  scale_y_continuous(limits = c(0, 40000),
                     labels = label_number(scale = 1e-3, suffix = "k")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        title = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.position = "right") + 
  coord_flip() +
  scale_fill_manual(values = c("black", "grey"))

fig1d
ggsave("fig1d.png", fig1d, width = 4.72, height = 2.47, units = "in", dpi = 300)

################################################################################
# Figure 1E: Transcriptional Complexity Across Tissues
################################################################################

# Prepare placenta data for multi-tissue comparison
classif$tissue <- "Placenta"
expression.placenta <- classif[, c("isoform", "associated_gene", "tissue", "FL")]
names(expression.placenta) <- c("transcript_id", "gene_id", "tissue", "expression")

# Combine GTEx and placenta data
expression_filtered <- rbind(expression_filtered, expression.placenta)
expression_filtered <- data.frame(expression_filtered)

# Calculate transcripts per gene for each tissue
tissues <- c(tissues, "Placenta")
txbyg_freq <- data.frame(matrix(ncol = 3, nrow = 0))
names(txbyg_freq) <- c("gene_id", "Freq", "tissue")

for(t in tissues) {
  df <- data.frame(table(expression_filtered[expression_filtered$tissue == t, "gene_id"]))
  df$tissue <- t
  names(df)[1] <- "gene_id"
  txbyg_freq <- rbind(txbyg_freq, df)
}

# Create box plot comparing transcriptional complexity across tissues
fig1e <- ggplot(txbyg_freq, aes(x = tissue, y = Freq)) +
  geom_boxplot(outliers = FALSE) +
  stat_summary(fun.y = mean, geom = "point", shape = 20, size = 4, 
               color = "red", fill = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        axis.title.y = element_blank(), 
        title = element_text(size = 8),
        axis.title.x = element_blank()) + 
  coord_flip()

fig1e
ggsave("fig1e.png", fig1e, width = 4.72, height = 2.47, units = "in", dpi = 300)

################################################################################
# End of Analysis
################################################################################