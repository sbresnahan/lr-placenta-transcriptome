################################################################################
# Short-Read Transcriptomics Analysis: Assembly vs Reference Comparison
# Author: Sean T. Bresnahan
# Description: Comprehensive comparison of short-read RNA-seq quantification using
#              the filtered placental transcriptome assembly versus GENCODE v45
#              reference. Analyzes expression patterns, statistical properties,
#              quantification uncertainty, transcript features, and coding potential
#              across two independent cohorts (GUSTO N=200, Gen3G N=152).
#
# Input Files:
#   - Filtered assembly classification (ESPRESSO_corrected_SC_filtered_classification.txt)
#   - Original assembly classification (ESPRESSO_corrected_classification.txt)
#   - Salmon quantification results (assembly and GENCODE for both cohorts)
#   - Sample metadata (GUSTO and Gen3G cohort information)
#   - Transcript-to-gene mappings (assembly and GENCODE GTF files)
#   - Coding potential analysis results (CPC2.txt)
#   - ORF predictions (longest_orfs.gff3)
#   - BLASTP results (blastp.outfmt6, ESPRESSO_corrected_SC_filtered_blastp.outfmt6)
#
# Output Figures:
#   - Main Manuscript: 1G, 1H, 1I, 2C
#     * 1G: Total expression comparison across cohorts and annotations
#     * 1H: Mean expression distribution by transcript categories (GUSTO)
#     * 1I: Inferential uncertainty by structural category
#     * 2C: Coding isoforms with CDS support from placenta MS/MS
#   - Supplementary: S4A-S4G, S5A-S5E, S6A-S6D
#     * S4A: PCA of GUSTO samples for quality control
#     * S4B: Mean/variance expression patterns (Gen3G)
#     * S4C-S4E: Inferential uncertainty analysis (both cohorts)
#     * S4F-S4G: Pairwise sample correlations by transcript category
#     * S5A-S5E: Transcript feature analysis (FL reads, length, exons, correlations)
#     * S6A-S6D: Coding potential and ORF analysis
################################################################################

################################################################################
# Environment Setup
################################################################################

# Load required libraries
library(dplyr)         # Data manipulation
library(rtracklayer)   # Genomic data import/export
library(ggplot2)       # Data visualization
library(tidyr)         # Data tidying
library(cowplot)       # Plot arrangement
library(plyr)          # Data manipulation (legacy functions)
library(scales)        # Scale functions for ggplot2
library(ggExtra)       # Additional ggplot2 functionality
library(edgeR)         # RNA-seq analysis toolkit
library(tximport)      # Transcript abundance import and summarization
library(DESeq2)        # Differential expression analysis
library(fishpond)      # Transcript-level inferential variance
library(tximeta)       # Metadata-aware transcript quantification import
library(ggdist)        # Distribution visualization
library(lme4)          # Mixed-effects models
library(lmerTest)      # Tests for mixed-effects models
library(foreach)       # Parallel processing support
library(doParallel)    # Parallel backend
library(emmeans)       # Estimated marginal means
library(enrichR)       # Gene enrichment analysis
library(org.Hs.eg.db)  # Human genome annotation

################################################################################
# Global Configuration
################################################################################

# Define color palettes for consistent visualization
myPalette <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#ffd92f",
               "#e5c494", "#c87774", "#d3adaf", "#b3b3b3")
myPalette2 <- c(myPalette[1:4], "grey20")

################################################################################
# Data Import and Processing
################################################################################

# Import SQANTI3 classification for filtered assembly
classif <- read.table("ESPRESSO_corrected_SC_filtered_classification.txt", header = TRUE)

# Standardize structural category definitions and factor levels
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

# Standardize column naming
names(classif)[1] <- "transcript_id"

# Import original assembly classification for comparison
classif.original <- read.table("/Users/stbresnahan/Desktop/Manuscripts/Placenta/ESPRESSO_original/ESPRESSO_corrected_classification.txt", header = TRUE)
classif.original$structural_category <- factor(as.character(classif.original$structural_category),
                                               levels = c("full-splice_match",
                                                          "incomplete-splice_match",
                                                          "novel_in_catalog",
                                                          "novel_not_in_catalog",
                                                          "genic",
                                                          "antisense",
                                                          "fusion",
                                                          "intergenic",
                                                          "genic_intron"))
classif.original$structural_category <- revalue(classif.original$structural_category, 
                                                c("full-splice_match" = "FSM",
                                                  "incomplete-splice_match" = "ISM",
                                                  "novel_in_catalog" = "NIC",
                                                  "novel_not_in_catalog" = "NNC",
                                                  "genic" = "Genic Genomic",
                                                  "antisense" = "Antisense",
                                                  "fusion" = "Fusion",
                                                  "intergenic" = "Intergenic",
                                                  "genic_intron" = "Genic Intron"))

# Use filtered classification as main dataset
classif.filtered <- classif

################################################################################
# GUSTO Cohort Data Processing
################################################################################

## Import and process GUSTO metadata
metadata <- read.csv("20220823-Full_200_RNAseq_covars_v3.csv")
metadata$SampleID <- paste0("J", metadata$ID)
metadata <- metadata[, c(51, 1:50)]
row.names(metadata) <- metadata$SampleID

# Add RNA quality metrics from separate file
metadata2 <- read.csv("20220823-Full_200_RNAseq_covars_v2.csv")
metadata$QC_RIN <- metadata2$QC_RIN
rm(metadata2)

# Prepare covariate data for statistical modeling
cordat <- metadata
cordat$GROUP <- factor(as.numeric(factor(cordat$GROUP)))
cordat$sex <- factor(as.numeric(factor(cordat$sex)))
cordat$gdm <- factor(cordat$gdm)
cordat$mother_ethnicity <- factor(cordat$mother_ethnicity)
cordat$htn_clean <- factor(as.numeric(factor(cordat$htn_clean)))
cordat$GA <- as.numeric(scale(cordat$GA))
cordat$c_weight_birth_kg <- scale(cordat$c_weight_birth_kg * 1000)

################################################################################
# Transcript-to-Gene Mapping for Filtered Assembly
################################################################################

# Import filtered assembly GTF and create transcript-to-gene mapping
tx2g.assembly <- data.frame(import("ESPRESSO_corrected_SC_filtered.gtf"))
tx2g.assembly <- tx2g.assembly[!duplicated(tx2g.assembly), ]
tx2g.assembly$score <- NULL
tx2g.assembly$phase <- NULL

# Handle novel transcripts without gene assignments
novelTx.assembly <- unique(tx2g.assembly[tx2g.assembly$gene_id == "NA", "transcript_id"])

# Use SQANTI3 associated_gene assignments for novel transcripts when available
novelTx.ascg <- classif[classif$transcript_id %in% novelTx.assembly, 
                        c("transcript_id", "associated_gene")]
novelTx.ascg <- novelTx.ascg[grep("ENSG", novelTx.ascg$associated_gene), ]

# Merge associated gene information
tx2g.assembly <- left_join(tx2g.assembly, novelTx.ascg, by = "transcript_id")
tx2g.assembly[is.na(tx2g.assembly$associated_gene), "associated_gene"] <- "NA"
tx2g.assembly[tx2g.assembly$gene_id == "NA", "gene_id"] <- 
  tx2g.assembly[tx2g.assembly$gene_id == "NA", "associated_gene"]
tx2g.assembly$associated_gene <- NULL

# Create gene IDs for remaining unannotated novel transcripts
novelTx.assembly.unannot <- setdiff(novelTx.assembly, novelTx.ascg$transcript_id)
tx2g.assembly[tx2g.assembly$transcript_id %in% novelTx.assembly.unannot, "gene_id"] <- 
  sapply(tx2g.assembly[tx2g.assembly$transcript_id %in% novelTx.assembly.unannot, "transcript_id"],
         function(x) paste(strsplit(x, ":")[[1]][1:3], collapse = ":"))

# Finalize transcript-to-gene mapping
tx2g.assembly <- tx2g.assembly[, c("transcript_id", "gene_id")]
tx2g.assembly <- tx2g.assembly[!duplicated(tx2g.assembly), ]

################################################################################
# GUSTO Assembly Quantification Processing
################################################################################

## Import Salmon quantification results and apply QU (quasi-UMI) correction
dirs.salmon.assembly <- list.dirs("COUNTS_assembly_SC", recursive = FALSE)

# Use fishpond's catchSalmon to identify overdispersion parameters
s.assembly <- catchSalmon(dirs.salmon.assembly)

# Create file list for tximport
files.list.assembly <- list()
for (i in 1:length(dirs.salmon.assembly)) {
  files.list.assembly[i] <- paste(unlist(dirs.salmon.assembly[i]), "quant.sf", sep = "/")
}
dirnames.assembly <- sapply(strsplit(unlist(dirs.salmon.assembly), "/"), "[",
                            length(unlist(strsplit(unlist(dirs.salmon.assembly), "/")[1])))
names(files.list.assembly) <- dirnames.assembly

# Import quantification data with QU correction for overdispersion
txi.assembly <- tximport(unlist(files.list.assembly), 
                         type = "salmon", 
                         tx2gene = tx2g.assembly,
                         geneIdCol = "gene_id", 
                         txIdCol = "transcript_id", 
                         dropInfReps = TRUE, 
                         txOut = TRUE)

# Apply QU correction to reduce technical variance
txi.assembly$counts <- txi.assembly$counts / s.assembly$annotation$Overdispersion
txi.assembly$abundance <- cpm(txi.assembly$counts)
txi.assembly.GUSTO <- txi.assembly

# Create DESeq2 dataset for downstream analysis
dds.assembly.GUSTO <- DESeqDataSetFromTximport(txi = txi.assembly.GUSTO, 
                                               colData = cordat, 
                                               design = ~1)

################################################################################
# Principal Component Analysis for Quality Control
################################################################################

# Variance stabilizing transformation for PCA
vst <- assay(varianceStabilizingTransformation(dds.assembly.GUSTO))

# Perform PCA to identify potential outliers and batch effects
pca <- prcomp(t(vst))
pca.df <- data.frame(pca$x[, 1:2])
pca.df$SampleID <- row.names(pca.df)
pca.df <- left_join(pca.df, cordat[, c("SampleID", "gdm")])
pca.df[is.na(pca.df)] <- 0

# Calculate variance explained by each PC
var_explained <- pca$sdev^2
pve <- var_explained / sum(var_explained)
pve[1:2] * 100  # PC1: 25.72%, PC2: 19.26%

## Figure S4A: PCA of GUSTO samples for outlier detection
figs4a <- ggplot(pca.df, aes(x = PC1, y = PC2, color = gdm)) +
  geom_point() +
  xlab("PC1: 25.72% variance") +
  ylab("PC2: 19.26% variance") +
  ggtitle("GUSTO Sample PCA")

figs4a
ggsave("figs4a.png", figs4a, width = 3, height = 3, units = "in", dpi = 300)

################################################################################
# GUSTO GENCODE Reference Quantification
################################################################################

# Import GENCODE v45 transcript-to-gene mapping
tx2g.gencode <- data.frame(import("gencode.v45.annotation.gtf"))
tx2g.gencode <- tx2g.gencode[tx2g.gencode$seqnames %in% 
                               c(paste0("chr", seq.int(1, 22, 1)), "chrX", "chrY", "chrM"), ]
tx2g.gencode <- tx2g.gencode[, c("transcript_id", "gene_id")]
tx2g.gencode <- tx2g.gencode[!is.na(tx2g.gencode$transcript_id), ]
tx2g.gencode <- tx2g.gencode[!duplicated(tx2g.gencode), ]

# Process GENCODE quantification data
dirs.salmon.gencode <- list.dirs("COUNTS_gencode", recursive = FALSE)
s.gencode <- catchSalmon(dirs.salmon.gencode)

files.list.gencode <- list()
for (i in 1:length(dirs.salmon.gencode)) {
  files.list.gencode[i] <- paste(unlist(dirs.salmon.gencode[i]), "quant.sf", sep = "/")
}
dirnames.gencode <- sapply(strsplit(unlist(dirs.salmon.gencode), "/"), "[",
                           length(unlist(strsplit(unlist(dirs.salmon.gencode), "/")[1])))
names(files.list.gencode) <- dirnames.gencode

# Import and process GENCODE quantification with QU correction
txi.gencode <- tximport(unlist(files.list.gencode), 
                        type = "salmon", 
                        tx2gene = tx2g.gencode,
                        geneIdCol = "gene_id", 
                        txIdCol = "transcript_id", 
                        dropInfReps = TRUE, 
                        txOut = TRUE)

txi.gencode$counts <- txi.gencode$counts / s.gencode$annotation$Overdispersion
txi.gencode$abundance <- cpm(txi.gencode$counts)
txi.gencode.GUSTO <- txi.gencode

dds.gencode.GUSTO <- DESeqDataSetFromTximport(txi = txi.gencode.GUSTO, 
                                              colData = cordat, 
                                              design = ~1)

################################################################################
# Gen3G Cohort Data Processing
################################################################################

## Import and process Gen3G metadata from multiple dbGaP files
metadata <- read.table(gzfile("phs003151.v1.pht013460.v1.p1.c1.Gen3G_Pooled_Visit2_Subject_Phenotypes.HMB-IRB.txt.gz"), 
                       header = TRUE)

metadata2 <- read.table(gzfile("phs003151.v1.pht013462.v1.p1.c1.Gen3G_RNASEQ_Sample_Attributes.HMB-IRB.txt.gz"), 
                        header = TRUE, sep = "\t")
metadata2 <- metadata2[, c(2, 1, 5:7)]

metadata3 <- read.csv("Gen3G_SraRunTable.csv")
metadata3 <- metadata3[, c("Run", "biospecimen_repository_sample_id", "submitted_subject_id")]
names(metadata3)[c(2, 3)] <- c("SAMPLE_ID", "SUBJECT_ID")

# Merge metadata files
metadata3 <- left_join(metadata3, metadata2, by = "SAMPLE_ID")
metadata3 <- metadata3[!is.na(metadata3$dbGaP_Sample_ID), ]
metadata <- left_join(metadata3, metadata, by = "SUBJECT_ID")

# Process subject-level sex information (handle multiple entries per subject)
metadata4 <- read.table(gzfile("phs003151.v1.pht013453.v1.p1.Gen3G_Subject.MULTI.txt.gz"), 
                        header = TRUE, sep = "\t")
metadata4 <- metadata4[, c(5, 4)]
names(metadata4)[1] <- "SUBJECT_ID"
metadata4 <- metadata4[!duplicated(metadata4), ]
metadata4 <- metadata4[-2, ]  # Remove header duplicate
metadata4 <- metadata4[!grepl(";", metadata4$SUBJECT_ID), ]

# Resolve multiple sex entries per subject (prioritize male if any male entry exists)
metadata5 <- data.frame(matrix(ncol = 2, nrow = 0))
names(metadata5) <- names(metadata4)
sids <- unique(metadata4$SUBJECT_ID)

for(i in 1:length(sids)) {
  x <- metadata4[metadata4$SUBJECT_ID == sids[i], ]
  if(any(x$SEX == "M")) {
    metadata5 <- rbind(metadata5, data.frame(SUBJECT_ID = sids[i], SEX = "M"))
  } else {
    metadata5 <- rbind(metadata5, data.frame(SUBJECT_ID = sids[i], SEX = "F"))
  }
}
rm(x)

metadata <- left_join(metadata, metadata5, by = "SUBJECT_ID")

# Filter for fetal-facing placental samples only
metadata <- metadata[metadata$HISTOLOGICAL_TYPE == "Fetal facing side", ]

# Add birth weight information
metadata6 <- read.table(gzfile("phs003151.v1.pht013457.v1.p1.c1.Gen3G_Pooled_Delivery_Subject_Phenotypes.HMB-IRB.txt.gz"),
                        header = TRUE, sep = "\t")[, c("SUBJECT_ID", "Accouchement_Poids")]
metadata <- left_join(metadata, metadata6, by = "SUBJECT_ID")

# Clean up temporary objects
rm(metadata2, metadata3, metadata4, metadata5, metadata6)

# Prepare Gen3G covariate data
cordat <- metadata
names(cordat)[which(names(cordat) == "Visite_Sem_Gestation_V2")] <- "GA"
cordat$GA <- as.numeric(scale(cordat$GA))

names(cordat)[which(names(cordat) == "SEX")] <- "sex"
cordat$sex <- factor(cordat$sex)

names(cordat)[which(names(cordat) == "Visite_Diagnostic_de_DG_V2")] <- "gdm"
cordat <- cordat[!cordat$gdm == 999, ]
cordat$gdm <- factor(cordat$gdm)

names(cordat)[31] <- "birth_weight"
cordat$birth_weight <- scale(cordat$birth_weight)

################################################################################
# Gen3G Quantification Processing (Assembly and GENCODE)
################################################################################

## Gen3G Assembly quantification
dirs.salmon.assembly <- list.dirs("COUNTS_assembly_SC", 
                                  recursive = FALSE)
names(dirs.salmon.assembly) <- basename(dirs.salmon.assembly)
dirs.salmon.assembly <- dirs.salmon.assembly[names(dirs.salmon.assembly) %in% cordat$Run]

s.assembly <- catchSalmon(dirs.salmon.assembly)

files.list.assembly <- list()
for (i in 1:length(dirs.salmon.assembly)) {
  files.list.assembly[i] <- paste(unlist(dirs.salmon.assembly[i]), "quant.sf", sep = "/")
}
dirnames.assembly <- sapply(strsplit(unlist(dirs.salmon.assembly), "/"), "[",
                            length(unlist(strsplit(unlist(dirs.salmon.assembly), "/")[1])))
names(files.list.assembly) <- dirnames.assembly

# Process Gen3G assembly quantification
txi.assembly <- tximport(unlist(files.list.assembly), 
                         type = "salmon", 
                         tx2gene = tx2g.assembly,
                         geneIdCol = "gene_id", 
                         txIdCol = "transcript_id", 
                         dropInfReps = TRUE, 
                         txOut = TRUE)

txi.assembly$counts <- txi.assembly$counts / s.assembly$annotation$Overdispersion
txi.assembly$abundance <- cpm(txi.assembly$counts)
txi.assembly.Gen3G <- txi.assembly
rm(txi.assembly)

dds.assembly.Gen3G <- DESeqDataSetFromTximport(txi = txi.assembly.Gen3G, 
                                               colData = cordat, 
                                               design = ~1)

## Gen3G GENCODE quantification
dirs.salmon.gencode <- list.dirs("COUNTS_gencode", 
                                 recursive = FALSE)
names(dirs.salmon.gencode) <- basename(dirs.salmon.gencode)
dirs.salmon.gencode <- dirs.salmon.gencode[names(dirs.salmon.gencode) %in% cordat$Run]

s.gencode <- catchSalmon(dirs.salmon.gencode)

files.list.gencode <- list()
for (i in 1:length(dirs.salmon.gencode)) {
  files.list.gencode[i] <- paste(unlist(dirs.salmon.gencode[i]), "quant.sf", sep = "/")
}
dirnames.gencode <- sapply(strsplit(unlist(dirs.salmon.gencode), "/"), "[",
                           length(unlist(strsplit(unlist(dirs.salmon.gencode), "/")[1])))
names(files.list.gencode) <- dirnames.gencode

txi.gencode <- tximport(unlist(files.list.gencode), 
                        type = "salmon", 
                        tx2gene = tx2g.gencode,
                        geneIdCol = "gene_id", 
                        txIdCol = "transcript_id", 
                        dropInfReps = TRUE, 
                        txOut = TRUE)

txi.gencode$counts <- txi.gencode$counts / s.gencode$annotation$Overdispersion
txi.gencode$abundance <- cpm(txi.gencode$counts)
txi.gencode.Gen3G <- txi.gencode
rm(txi.gencode)

dds.gencode.Gen3G <- DESeqDataSetFromTximport(txi = txi.gencode.Gen3G, 
                                              colData = cordat, 
                                              design = ~1)

################################################################################
# Figure 1G: Total Expression Comparison Across Cohorts and Annotations
################################################################################

## Calculate normalized TPM-like expression values for comparison
# GUSTO Assembly
fig1g.dat.1 <- counts(estimateSizeFactors(dds.assembly.GUSTO), normalized = TRUE)
fig1g.dat.1 <- fig1g.dat.1 / txi.assembly.GUSTO$length
fig1g.dat.1 <- log(t(t(fig1g.dat.1) * 1e6 / colSums(fig1g.dat.1)) + 1)
fig1g.dat.1 <- data.frame(colSums(fig1g.dat.1))
fig1g.dat.1$sample_id <- row.names(fig1g.dat.1)
fig1g.dat.1 <- fig1g.dat.1[, c(2, 1)]
names(fig1g.dat.1)[2] <- "assembly"

# GUSTO GENCODE
fig1g.dat.2 <- counts(estimateSizeFactors(dds.gencode.GUSTO), normalized = TRUE)
fig1g.dat.2 <- fig1g.dat.2 / txi.gencode.GUSTO$length
fig1g.dat.2 <- log(t(t(fig1g.dat.2) * 1e6 / colSums(fig1g.dat.2)) + 1)
fig1g.dat.2 <- data.frame(colSums(fig1g.dat.2))
fig1g.dat.2$sample_id <- row.names(fig1g.dat.2)
fig1g.dat.2 <- fig1g.dat.2[, c(2, 1)]
names(fig1g.dat.2)[2] <- "gencode"

# Combine GUSTO data
fig1g.dat.a <- data.frame(sample_id = fig1g.dat.1$sample_id)
fig1g.dat.a <- left_join(fig1g.dat.a, fig1g.dat.1)
fig1g.dat.a <- left_join(fig1g.dat.a, fig1g.dat.2)

# Gen3G Assembly
fig1g.dat.3 <- counts(estimateSizeFactors(dds.assembly.Gen3G), normalized = TRUE)
fig1g.dat.3 <- fig1g.dat.3 / txi.assembly.Gen3G$length
fig1g.dat.3 <- log(t(t(fig1g.dat.3) * 1e6 / colSums(fig1g.dat.3)) + 1)
fig1g.dat.3 <- data.frame(colSums(fig1g.dat.3))
fig1g.dat.3$sample_id <- row.names(fig1g.dat.3)
fig1g.dat.3 <- fig1g.dat.3[, c(2, 1)]
names(fig1g.dat.3)[2] <- "assembly"

# Gen3G GENCODE
fig1g.dat.4 <- counts(estimateSizeFactors(dds.gencode.Gen3G), normalized = TRUE)
fig1g.dat.4 <- fig1g.dat.4 / txi.gencode.Gen3G$length
fig1g.dat.4 <- log(t(t(fig1g.dat.4) * 1e6 / colSums(fig1g.dat.4)) + 1)
fig1g.dat.4 <- data.frame(colSums(fig1g.dat.4))
fig1g.dat.4$sample_id <- row.names(fig1g.dat.4)
fig1g.dat.4 <- fig1g.dat.4[, c(2, 1)]
names(fig1g.dat.4)[2] <- "gencode"

# Combine Gen3G data
fig1g.dat.b <- data.frame(sample_id = fig1g.dat.3$sample_id)
fig1g.dat.b <- left_join(fig1g.dat.b, fig1g.dat.3)
fig1g.dat.b <- left_join(fig1g.dat.b, fig1g.dat.4)

# Create comprehensive comparison dataset
fig1g.dat <- data.frame(
  Cohort = c(rep.int("GUSTO n=200", length(fig1g.dat.a$assembly) * 2),
             rep.int("Gen3G n=152", length(fig1g.dat.b$assembly) * 2)),
  annotation = c(rep.int("lr-assembly", length(fig1g.dat.a$assembly)),
                 rep.int("GENCODEv45", length(fig1g.dat.a$gencode)),
                 rep.int("lr-assembly", length(fig1g.dat.b$assembly)),
                 rep.int("GENCODEv45", length(fig1g.dat.b$gencode))),
  sum = c(fig1g.dat.a$assembly, fig1g.dat.a$gencode, 
          fig1g.dat.b$assembly, fig1g.dat.b$gencode)
)

## Generate visualization
fig1g <- ggplot(fig1g.dat, aes(x = Cohort, y = sum, fill = annotation)) +
  geom_violin(alpha = 0.5, position = position_dodge(0.9)) + 
  geom_boxplot(alpha = 0.5, width = 0.1, position = position_dodge(0.9)) + 
  theme_minimal() +
  theme(legend.position = "top",
        axis.title.x = element_blank(),
        legend.title = element_blank()) +
  labs(y = expression("TPM x10³")) +
  ggtitle("Total expression across isoforms per sample") +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values = c("royalblue", "seagreen"))

fig1g
ggsave("fig1g.png", fig1g, width = 4.5, height = 4, units = "in", dpi = 300)

## Statistical analysis of expression differences
fig1g.dat$sample <- c(paste0("GUSTO_", 1:200), paste0("GUSTO_", 1:200),
                      paste0("Gen3G_", 1:152), paste0("Gen3G_", 1:152))

# Mixed-effects model accounting for paired samples
model <- lmer(sum ~ annotation + Cohort + (1 | sample), data = fig1g.dat)
summary(model)

################################################################################
# Expression Pattern Analysis: Mean and Variance by Transcript Categories
################################################################################

## Define transcript categories for comparison
txs.all <- union(tx2g.assembly$transcript_id, tx2g.gencode$transcript_id)
txs.both <- intersect(tx2g.assembly$transcript_id, tx2g.gencode$transcript_id)

# Load raw assembly data to identify filtering effects
gtf.assembly_all <- "ESPRESSO.gtf"
all_assembled_tx <- unique(data.frame(import(gtf.assembly_all))$transcript_id)

# Classify transcripts by assembly/filtering status
not_assembled <- setdiff(txs.all, all_assembled_tx)
filtered <- setdiff(all_assembled_tx, tx2g.assembly$transcript_id)
shared <- txs.both
novel <- setdiff(tx2g.assembly$transcript_id, txs.both)

################################################################################
# Figure 1H: Mean Expression Distribution by Transcript Categories (GUSTO)
################################################################################

## Calculate mean expression for each transcript across GUSTO samples
fig1h1.dat.1 <- counts(estimateSizeFactors(dds.assembly.GUSTO), normalized = TRUE)
fig1h1.dat.1 <- data.frame(rowMeans(fig1h1.dat.1))
fig1h1.dat.1$transcript_id <- row.names(fig1h1.dat.1)
fig1h1.dat.1 <- fig1h1.dat.1[, c(2, 1)]
names(fig1h1.dat.1)[2] <- "lr-assembly"

fig1h1.dat.2 <- counts(estimateSizeFactors(dds.gencode.GUSTO), normalized = TRUE)
fig1h1.dat.2 <- data.frame(rowMeans(fig1h1.dat.2))
fig1h1.dat.2$transcript_id <- row.names(fig1h1.dat.2)
fig1h1.dat.2 <- fig1h1.dat.2[, c(2, 1)]
names(fig1h1.dat.2)[2] <- "GENCODEv45"

# Combine data and classify transcripts
fig1h1.dat.a <- data.frame(transcript_id = txs.all)
fig1h1.dat.a <- left_join(fig1h1.dat.a, fig1h1.dat.1)
fig1h1.dat.a <- left_join(fig1h1.dat.a, fig1h1.dat.2)
fig1h1.dat.a[is.na(fig1h1.dat.a)] <- 0

# Assign category labels
fig1h1.dat.a$set <- NA
fig1h1.dat.a[fig1h1.dat.a$transcript_id %in% shared, "set"] <- "shared"
fig1h1.dat.a[fig1h1.dat.a$transcript_id %in% novel, "set"] <- "novel"
fig1h1.dat.a[fig1h1.dat.a$transcript_id %in% filtered, "set"] <- "filtered"
fig1h1.dat.a[fig1h1.dat.a$transcript_id %in% not_assembled, "set"] <- "absent"

# Convert to long format for visualization
fig1h1.dat.a <- fig1h1.dat.a[, -1] %>%
  pivot_longer(
    cols = c(`lr-assembly`, GENCODEv45),
    names_to = "annotation",
    values_to = "mean"
  )

# Remove impossible combinations (novel transcripts can't be in GENCODE, etc.)
fig1h1.dat.a$key <- paste0(fig1h1.dat.a$set, fig1h1.dat.a$annotation)
fig1h1.dat.a <- fig1h1.dat.a[!fig1h1.dat.a$key == "novelGENCODEv45", ]
fig1h1.dat.a <- fig1h1.dat.a[!fig1h1.dat.a$key == "absentlr-assembly", ]
fig1h1.dat.a <- fig1h1.dat.a[!fig1h1.dat.a$key == "filteredlr-assembly", ]

fig1h1.dat.a$set <- factor(fig1h1.dat.a$set, levels = c("shared", "novel", "filtered", "absent"))
fig1h1.dat.a <- fig1h1.dat.a[order(fig1h1.dat.a$set), ]
fig1h1.dat.a$key <- NULL
fig1h1.dat.a <- fig1h1.dat.a[complete.cases(fig1h1.dat.a), ]

## Create visualization
fig1h1.a <- ggplot(fig1h1.dat.a,
                   aes(x = set, y = mean + 1, fill = annotation)) +
  geom_boxplot(alpha = 0.5, outliers = FALSE) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(),
        panel.grid.minor.y = element_line(),
        legend.title = element_blank(),
        legend.position = "top") +
  xlab("Isoforms in lr-assembly vs GENCODEv45") +
  ylab("Mean") +
  scale_fill_manual(values = c("royalblue", "seagreen")) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  ggtitle("Mean of short-read counts (GUSTO n=200)")

fig1h1.a
ggsave("fig1h1.png", fig1h1.a, width = 5, height = 5, units = "in", dpi = 300)

################################################################################
# Figure S4B: Mean Expression Distribution (Gen3G)
################################################################################

## Repeat analysis for Gen3G cohort
figs4ba.dat.3 <- counts(estimateSizeFactors(dds.assembly.Gen3G), normalized = TRUE)
figs4ba.dat.3 <- data.frame(rowMeans(figs4ba.dat.3))
figs4ba.dat.3$transcript_id <- row.names(figs4ba.dat.3)
figs4ba.dat.3 <- figs4ba.dat.3[, c(2, 1)]
names(figs4ba.dat.3)[2] <- "lr-assembly"

figs4ba.dat.4 <- counts(estimateSizeFactors(dds.gencode.Gen3G), normalized = TRUE)
figs4ba.dat.4 <- data.frame(rowMeans(figs4ba.dat.4))
figs4ba.dat.4$transcript_id <- row.names(figs4ba.dat.4)
figs4ba.dat.4 <- figs4ba.dat.4[, c(2, 1)]
names(figs4ba.dat.4)[2] <- "GENCODEv45"

# Process Gen3G data with same classification scheme
figs4ba.dat.b <- data.frame(transcript_id = txs.all)
figs4ba.dat.b <- left_join(figs4ba.dat.b, figs4ba.dat.3)
figs4ba.dat.b <- left_join(figs4ba.dat.b, figs4ba.dat.4)
figs4ba.dat.b[is.na(figs4ba.dat.b)] <- 0

figs4ba.dat.b$set <- NA
figs4ba.dat.b[figs4ba.dat.b$transcript_id %in% shared, "set"] <- "shared"
figs4ba.dat.b[figs4ba.dat.b$transcript_id %in% novel, "set"] <- "novel"
figs4ba.dat.b[figs4ba.dat.b$transcript_id %in% filtered, "set"] <- "filtered"
figs4ba.dat.b[figs4ba.dat.b$transcript_id %in% not_assembled, "set"] <- "absent"

figs4ba.dat.b <- figs4ba.dat.b[, -1] %>%
  pivot_longer(
    cols = c(`lr-assembly`, GENCODEv45),
    names_to = "annotation",
    values_to = "mean"
  )

figs4ba.dat.b$key <- paste0(figs4ba.dat.b$set, figs4ba.dat.b$annotation)
figs4ba.dat.b <- figs4ba.dat.b[!figs4ba.dat.b$key == "novelGENCODEv45", ]
figs4ba.dat.b <- figs4ba.dat.b[!figs4ba.dat.b$key == "absentlr-assembly", ]
figs4ba.dat.b <- figs4ba.dat.b[!figs4ba.dat.b$key == "filteredlr-assembly", ]

figs4ba.dat.b$set <- factor(figs4ba.dat.b$set, levels = c("shared", "novel", "filtered", "absent"))
figs4ba.dat.b <- figs4ba.dat.b[order(figs4ba.dat.b$set), ]
figs4ba.dat.b$key <- NULL
figs4ba.dat.b <- figs4ba.dat.b[complete.cases(figs4ba.dat.b), ]

figs4ba.b <- ggplot(figs4ba.dat.b,
                    aes(x = set, y = mean + 1, fill = annotation)) +
  geom_boxplot(alpha = 0.5, outliers = FALSE) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(),
        panel.grid.minor.y = element_line(),
        legend.title = element_blank(),
        legend.position = "top") +
  xlab("Isoforms in lr-assembly vs GENCODEv45") +
  ylab("Mean") +
  scale_fill_manual(values = c("royalblue", "seagreen")) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  ggtitle("Mean of short-read counts (Gen3G n=152)")

figs4ba.b
ggsave("figs4ba.png", figs4ba.b, width = 5, height = 5, units = "in", dpi = 300)

################################################################################
# Statistical Analysis of Mean Expression Differences
################################################################################

## Combine both cohorts for statistical testing
fig1h1.dat.a <- data.frame(fig1h1.dat.a)
fig1h1.dat.a$cohort <- "GUSTO"
figs4ba.dat.b <- data.frame(figs4ba.dat.b)
figs4ba.dat.b$cohort <- "Gen3G"
fig1h1.dat.all <- rbind(fig1h1.dat.a, figs4ba.dat.b)

# Prepare factors for modeling
fig1h1.dat.all$set <- factor(fig1h1.dat.all$set, levels = c("shared", "novel", "filtered", "absent"))
fig1h1.dat.all$annotation <- factor(fig1h1.dat.all$annotation)
fig1h1.dat.all$cohort <- factor(fig1h1.dat.all$cohort)

# Nested linear model to test differences between categories
fit_nested <- lm(log1p(mean) ~ cohort:set + cohort:annotation, data = fig1h1.dat.all)
anova(fit_nested)

# Post-hoc pairwise comparisons
emm_nested <- emmeans(fit_nested, ~ set | cohort)
pairwise_contrasts <- contrast(emm_nested, method = "pairwise", adjust = "tukey")
summary(pairwise_contrasts, infer = c(TRUE, TRUE))

################################################################################
# Expression Variance Analysis
################################################################################

## Figure 1H: Variance Distribution by Transcript Categories (GUSTO)
fig1h2.dat.1 <- counts(estimateSizeFactors(dds.assembly.GUSTO), normalized = TRUE)
fig1h2.dat.1 <- data.frame(rowVars(fig1h2.dat.1))
fig1h2.dat.1$transcript_id <- row.names(fig1h2.dat.1)
fig1h2.dat.1 <- fig1h2.dat.1[, c(2, 1)]
names(fig1h2.dat.1)[2] <- "lr-assembly"

fig1h2.dat.2 <- counts(estimateSizeFactors(dds.gencode.GUSTO), normalized = TRUE)
fig1h2.dat.2 <- data.frame(rowVars(fig1h2.dat.2))
fig1h2.dat.2$transcript_id <- row.names(fig1h2.dat.2)
fig1h2.dat.2 <- fig1h2.dat.2[, c(2, 1)]
names(fig1h2.dat.2)[2] <- "GENCODEv45"

# Process variance data with same classification approach
fig1h2.dat.a <- data.frame(transcript_id = txs.all)
fig1h2.dat.a <- left_join(fig1h2.dat.a, fig1h2.dat.1)
fig1h2.dat.a <- left_join(fig1h2.dat.a, fig1h2.dat.2)
fig1h2.dat.a[is.na(fig1h2.dat.a)] <- 0

fig1h2.dat.a$set <- NA
fig1h2.dat.a[fig1h2.dat.a$transcript_id %in% shared, "set"] <- "shared"
fig1h2.dat.a[fig1h2.dat.a$transcript_id %in% novel, "set"] <- "novel"
fig1h2.dat.a[fig1h2.dat.a$transcript_id %in% filtered, "set"] <- "filtered"
fig1h2.dat.a[fig1h2.dat.a$transcript_id %in% not_assembled, "set"] <- "absent"

fig1h2.dat.a <- fig1h2.dat.a[, -1] %>%
  pivot_longer(
    cols = c(`lr-assembly`, GENCODEv45),
    names_to = "annotation",
    values_to = "mean"
  )

fig1h2.dat.a$key <- paste0(fig1h2.dat.a$set, fig1h2.dat.a$annotation)
fig1h2.dat.a <- fig1h2.dat.a[!fig1h2.dat.a$key == "novelGENCODEv45", ]
fig1h2.dat.a <- fig1h2.dat.a[!fig1h2.dat.a$key == "absentlr-assembly", ]
fig1h2.dat.a <- fig1h2.dat.a[!fig1h2.dat.a$key == "filteredlr-assembly", ]

fig1h2.dat.a$set <- factor(fig1h2.dat.a$set, levels = c("shared", "novel", "filtered", "absent"))
fig1h2.dat.a <- fig1h2.dat.a[order(fig1h2.dat.a$set), ]
fig1h2.dat.a$key <- NULL
fig1h2.dat.a <- fig1h2.dat.a[complete.cases(fig1h2.dat.a), ]

fig1h2.a <- ggplot(fig1h2.dat.a,
                   aes(x = set, y = mean + 1, fill = annotation)) +
  geom_boxplot(alpha = 0.5, outliers = FALSE) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(),
        panel.grid.minor.y = element_line(),
        legend.title = element_blank(),
        legend.position = "top") +
  xlab("Isoforms in lr-assembly vs GENCODEv45") +
  ylab("Variance") +
  scale_fill_manual(values = c("royalblue", "seagreen")) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  ggtitle("Variance of short-read counts (GUSTO n=200)")

fig1h2.a
ggsave("fig1h2.png", fig1h2.a, width = 5, height = 5, units = "in", dpi = 300)

## Figure S4B: Variance Distribution (Gen3G)
figs4bb.dat.3 <- counts(estimateSizeFactors(dds.assembly.Gen3G), normalized = TRUE)
figs4bb.dat.3 <- data.frame(rowVars(figs4bb.dat.3))
figs4bb.dat.3$transcript_id <- row.names(figs4bb.dat.3)
figs4bb.dat.3 <- figs4bb.dat.3[, c(2, 1)]
names(figs4bb.dat.3)[2] <- "lr-assembly"

figs4bb.dat.4 <- counts(estimateSizeFactors(dds.gencode.Gen3G), normalized = TRUE)
figs4bb.dat.4 <- data.frame(rowVars(figs4bb.dat.4))
figs4bb.dat.4$transcript_id <- row.names(figs4bb.dat.4)
figs4bb.dat.4 <- figs4bb.dat.4[, c(2, 1)]
names(figs4bb.dat.4)[2] <- "GENCODEv45"

figs4bb.dat.b <- data.frame(transcript_id = txs.all)
figs4bb.dat.b <- left_join(figs4bb.dat.b, figs4bb.dat.3)
figs4bb.dat.b <- left_join(figs4bb.dat.b, figs4bb.dat.4)
figs4bb.dat.b[is.na(figs4bb.dat.b)] <- 0

figs4bb.dat.b$set <- NA
figs4bb.dat.b[figs4bb.dat.b$transcript_id %in% shared, "set"] <- "shared"
figs4bb.dat.b[figs4bb.dat.b$transcript_id %in% novel, "set"] <- "novel"
figs4bb.dat.b[figs4bb.dat.b$transcript_id %in% filtered, "set"] <- "filtered"
figs4bb.dat.b[figs4bb.dat.b$transcript_id %in% not_assembled, "set"] <- "absent"

figs4bb.dat.b <- figs4bb.dat.b[, -1] %>%
  pivot_longer(
    cols = c(`lr-assembly`, GENCODEv45),
    names_to = "annotation",
    values_to = "mean"
  )

figs4bb.dat.b$key <- paste0(figs4bb.dat.b$set, figs4bb.dat.b$annotation)
figs4bb.dat.b <- figs4bb.dat.b[!figs4bb.dat.b$key == "novelGENCODEv45", ]
figs4bb.dat.b <- figs4bb.dat.b[!figs4bb.dat.b$key == "absentlr-assembly", ]
figs4bb.dat.b <- figs4bb.dat.b[!figs4bb.dat.b$key == "filteredlr-assembly", ]

figs4bb.dat.b$set <- factor(figs4bb.dat.b$set, levels = c("shared", "novel", "filtered", "absent"))
figs4bb.dat.b <- figs4bb.dat.b[order(figs4bb.dat.b$set), ]
figs4bb.dat.b$key <- NULL
figs4bb.dat.b <- figs4bb.dat.b[complete.cases(figs4bb.dat.b), ]

figs4bb.b <- ggplot(figs4bb.dat.b,
                    aes(x = set, y = mean + 1, fill = annotation)) +
  geom_boxplot(alpha = 0.5, outliers = FALSE) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(),
        panel.grid.minor.y = element_line(),
        legend.title = element_blank(),
        legend.position = "top") +
  xlab("Isoforms in lr-assembly vs GENCODEv45") +
  ylab("Variance") +
  scale_fill_manual(values = c("royalblue", "seagreen")) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  ggtitle("Variance of short-read counts (Gen3G n=152)")

figs4bb.b
ggsave("figs4bb.png", figs4bb.b, width = 5, height = 5, units = "in", dpi = 300)

## Statistical analysis of variance patterns
fig1h2 <- data.frame(fig1h2.dat.a)
fig1h2$cohort <- "GUSTO"
figs4bb.dat.b <- data.frame(figs4bb.dat.b)
figs4bb.dat.b$cohort <- "Gen3G"
fig1h2all <- rbind(fig1h2, figs4bb.dat.b)

fig1h2all$set <- factor(fig1h2all$set, levels = c("shared", "novel", "filtered", "absent"))
fig1h2all$annotation <- factor(fig1h2all$annotation)
fig1h2all$cohort <- factor(fig1h2all$cohort)
names(fig1h2all)[3] <- "variance"

fit_nested <- lm(log1p(variance) ~ cohort:set + cohort:annotation, data = fig1h2all)
anova(fit_nested)

emm_nested <- emmeans(fit_nested, ~ set | cohort)
pairwise_contrasts <- contrast(emm_nested, method = "pairwise", adjust = "tukey")
summary(pairwise_contrasts, infer = c(TRUE, TRUE))

################################################################################
# Inferential Uncertainty Analysis
################################################################################

## Prepare classification data for inferential variance analysis
classif.pcorr <- classif
classif.pcorr$structural_category <- revalue(classif.pcorr$structural_category, 
                                             c("Genic Genomic" = "Other",
                                               "Antisense" = "Other",
                                               "Fusion" = "Other",
                                               "Intergenic" = "Other",
                                               "Genic Intron" = "Other"))
names(classif.pcorr)[1] <- "isoform"

################################################################################
# GUSTO Inferential Uncertainty Analysis
################################################################################

## Setup data paths and metadata for tximeta
DIR_COUNTS_assembly <- "COUNTS_assembly_SC"
DIR_COUNTS_gencode <- "COUNTS_gencode"

# Reload GUSTO metadata for inferential variance analysis
metadata <- read.csv("20220823-Full_200_RNAseq_covars_v3.csv")
metadata$SampleID <- paste0("J", metadata$ID)
metadata <- metadata[, c(51, 1:50)]
row.names(metadata) <- metadata$SampleID

metadata2 <- read.csv("20220823-Full_200_RNAseq_covars_v2.csv")
metadata$QC_RIN <- metadata2$QC_RIN
rm(metadata2)

cordat <- metadata
cordat$GROUP <- factor(as.numeric(factor(cordat$GROUP)))
cordat$sex <- factor(as.numeric(factor(cordat$sex)))
cordat$gdm <- factor(cordat$gdm)
cordat$mother_ethnicity <- factor(cordat$mother_ethnicity)
cordat$htn_clean <- factor(as.numeric(factor(cordat$htn_clean)))
cordat$GA <- as.numeric(scale(cordat$GA))
cordat$c_weight_birth_kg <- scale(cordat$c_weight_birth_kg * 1000)

# Prepare file paths for tximeta
coldata.x <- cordat
coldata.x$files <- paste0(DIR_COUNTS_assembly, "/", coldata.x$SampleID, "/quant.sf")
names(coldata.x)[1] <- "names"

coldata.y <- cordat
coldata.y$files <- paste0(DIR_COUNTS_gencode, "/", coldata.y$SampleID, "/quant.sf")
names(coldata.y)[1] <- "names"

## Assembly inferential variance calculation
setTximetaBFC("ESPRESSO_corrected_SC_filtered")
makeLinkedTxome(indexDir = "ESPRESSO_corrected_SC_filtered",
                source = "ESPRESSO assembly of placenta lr-rnaseq (filtered)",
                organism = "Homo sapiens",
                genome = "GRCh38",
                release = "p14",
                fasta = "ESPRESSO_corrected_SC_filtered/salmon_index/gentrome.fa",
                gtf = "ESPRESSO_corrected_SC_filtered.gtf",
                write = FALSE)

# Load data with inferential replicates
se.x <- tximeta(coldata.x[, c(53, 1:52)], type = "salmon")
se.x <- labelKeep(se.x)
se.x <- se.x[mcols(se.x)$keep, ]

# Calculate inferential relative variance (InfRV)
infrv.assembly.GUSTO <- computeInfRV(se.x)
infrv.assembly.GUSTO <- mcols(infrv.assembly.GUSTO)$meanInfRV
rm(se.x)
gc()

## GENCODE inferential variance calculation
setTximetaBFC("GENCODE")
makeLinkedTxome(indexDir = "GENCODE/gencode.v45.salmon_index/gencode_v45",
                source = "GENCODEv45",
                organism = "Homo sapiens",
                genome = "GRCh38",
                release = "p14",
                fasta = "GENCODE/gencode.v45.salmon_index/gentrome.fa",
                gtf = "gencode_v45_annotation.gtf",
                write = FALSE)

se.y <- tximeta(coldata.y[, c(53, 1:52)], type = "salmon")
se.y <- labelKeep(se.y)
se.y <- se.y[mcols(se.y)$keep, ]

infrv.gencode.GUSTO <- computeInfRV(se.y)
infrv.gencode.GUSTO <- mcols(infrv.gencode.GUSTO)$meanInfRV
rm(se.y)
gc()

## Process and visualize GUSTO inferential uncertainty
infrv.GUSTO <- data.frame(transcript_id = union(names(infrv.assembly.GUSTO),
                                                names(infrv.gencode.GUSTO)))

infrv.assembly.GUSTO <- data.frame(infrv.assembly.GUSTO)
infrv.assembly.GUSTO$transcript_id <- row.names(infrv.assembly.GUSTO)
infrv.assembly.GUSTO <- infrv.assembly.GUSTO[, c(2, 1)]
names(infrv.assembly.GUSTO)[2] <- "lr-assembly"
infrv.GUSTO <- left_join(infrv.GUSTO, infrv.assembly.GUSTO)

infrv.gencode.GUSTO <- data.frame(infrv.gencode.GUSTO)
infrv.gencode.GUSTO$transcript_id <- row.names(infrv.gencode.GUSTO)
infrv.gencode.GUSTO <- infrv.gencode.GUSTO[, c(2, 1)]
names(infrv.gencode.GUSTO)[2] <- "GENCODEv45"
infrv.GUSTO <- left_join(infrv.GUSTO, infrv.gencode.GUSTO)

# Convert to long format and add structural categories
infrv.GUSTO <- infrv.GUSTO %>%
  pivot_longer(
    cols = c(`lr-assembly`, GENCODEv45),
    names_to = "annotation",
    values_to = "InfRV"
  )

infrv.GUSTO <- infrv.GUSTO[complete.cases(infrv.GUSTO), ]
names(infrv.GUSTO)[1] <- "isoform"

infrv.GUSTO <- left_join(infrv.GUSTO, classif.pcorr[, c("isoform", "structural_category")])
infrv.GUSTO$structural_category <- as.character(infrv.GUSTO$structural_category)
infrv.GUSTO[infrv.GUSTO$annotation == "GENCODEv45", "structural_category"] <- "GENCODE"
infrv.GUSTO$structural_category <- factor(infrv.GUSTO$structural_category,
                                          levels = c("GENCODE", "FSM", "ISM", "NIC", "NNC", "Other"))

## Figure S4C: Inferential uncertainty distribution (GUSTO)
figs4c <- ggplot(infrv.GUSTO, 
                 aes(x = log10(InfRV), fill = annotation, y = after_stat(density))) +
  geom_histogram(color = "#e9ecef", alpha = 0.4, position = 'identity', bins = 60) +
  scale_fill_manual(values = c("royalblue", "seagreen")) +
  theme_classic() +
  labs(fill = "") +
  ggtitle("") +
  geom_vline(xintercept = mean(log10(as.numeric(unlist(infrv.GUSTO[infrv.GUSTO$annotation == "lr-assembly", "InfRV"])))),
             color = "seagreen", linewidth = 1, linetype = "dashed") +
  geom_vline(xintercept = mean(log10(as.numeric(unlist(infrv.GUSTO[infrv.GUSTO$annotation == "GENCODEv45", "InfRV"])))),
             color = "royalblue", linewidth = 1, linetype = "dashed") +
  theme(panel.grid.major.y = element_line(),
        panel.grid.minor.y = element_line(),
        panel.grid.major.x = element_line(),
        panel.grid.minor.x = element_line(),
        legend.position = "top") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("mean InfRV") + 
  ylab("Density") +
  ggtitle("Isoform-level inferential uncertainty\nAll isoforms (GUSTO n=200)")

figs4c
ggsave("figs4c.png", figs4c, width = 6, height = 6, units = "in", dpi = 300)

## Figure 1I: Inferential uncertainty by structural category (GUSTO)
fig1i <- ggplot(infrv.GUSTO, aes(x = structural_category, fill = structural_category, y = scale(log10(InfRV)))) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.5, width = 0.1) +
  stat_halfeye(adjust = .5, width = .6, justification = -.2, .width = 0, point_colour = NA) +
  theme_minimal() + 
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.title.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_fill_manual(values = c("goldenrod", myPalette2)) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("mean InfRV") +
  ggtitle("Isoform-level inferential uncertainty")

fig1i
ggsave("fig1i.png", fig1i, width = 4, height = 3, units = "in", dpi = 300)

################################################################################
# Gen3G Inferential Uncertainty Analysis
################################################################################

## Setup Gen3G metadata for inferential variance
DIR_COUNTS_assembly <- "COUNTS_assembly_SC"
DIR_COUNTS_gencode <- "COUNTS_gencode"

# Reload and process Gen3G metadata (same as before)
metadata <- read.table(gzfile("phs003151.v1.pht013460.v1.p1.c1.Gen3G_Pooled_Visit2_Subject_Phenotypes.HMB-IRB.txt.gz"), 
                       header = TRUE)

metadata2 <- read.table(gzfile("phs003151.v1.pht013462.v1.p1.c1.Gen3G_RNASEQ_Sample_Attributes.HMB-IRB.txt.gz"), 
                        header = TRUE, sep = "\t")
metadata2 <- metadata2[, c(2, 1, 5:7)]

metadata3 <- read.csv("Gen3G_SraRunTable.csv")
metadata3 <- metadata3[, c("Run", "biospecimen_repository_sample_id", "submitted_subject_id")]
names(metadata3)[c(2, 3)] <- c("SAMPLE_ID", "SUBJECT_ID")

metadata3 <- left_join(metadata3, metadata2, by = "SAMPLE_ID")
metadata3 <- metadata3[!is.na(metadata3$dbGaP_Sample_ID), ]
metadata <- left_join(metadata3, metadata, by = "SUBJECT_ID")

metadata4 <- read.table(gzfile("phs003151.v1.pht013453.v1.p1.Gen3G_Subject.MULTI.txt.gz"), 
                        header = TRUE, sep = "\t")
metadata4 <- metadata4[, c(5, 4)]
names(metadata4)[1] <- "SUBJECT_ID"
metadata4 <- metadata4[!duplicated(metadata4), ]
metadata4 <- metadata4[-2, ]
metadata4 <- metadata4[!grepl(";", metadata4$SUBJECT_ID), ]

metadata5 <- data.frame(matrix(ncol = 2, nrow = 0))
names(metadata5) <- names(metadata4)
sids <- unique(metadata4$SUBJECT_ID)

for(i in 1:length(sids)) {
  x <- metadata4[metadata4$SUBJECT_ID == sids[i], ]
  if(any(x$SEX == "M")) {
    metadata5 <- rbind(metadata5, data.frame(SUBJECT_ID = sids[i], SEX = "M"))
  } else {
    metadata5 <- rbind(metadata5, data.frame(SUBJECT_ID = sids[i], SEX = "F"))
  }
}
rm(x, sids)

metadata <- left_join(metadata, metadata5, by = "SUBJECT_ID")
metadata <- metadata[metadata$HISTOLOGICAL_TYPE == "Fetal facing side", ]

metadata6 <- read.table(gzfile("phs003151.v1.pht013457.v1.p1.c1.Gen3G_Pooled_Delivery_Subject_Phenotypes.HMB-IRB.txt.gz"),
                        header = TRUE, sep = "\t")[, c("SUBJECT_ID", "Accouchement_Poids")]
metadata <- left_join(metadata, metadata6, by = "SUBJECT_ID")

rm(metadata2, metadata3, metadata4, metadata5, metadata6)

cordat <- metadata
names(cordat)[which(names(cordat) == "Visite_Sem_Gestation_V2")] <- "GA"
cordat$GA <- as.numeric(scale(cordat$GA))

names(cordat)[which(names(cordat) == "SEX")] <- "sex"
cordat$sex <- factor(cordat$sex)

names(cordat)[which(names(cordat) == "Visite_Diagnostic_de_DG_V2")] <- "gdm"
cordat <- cordat[!cordat$gdm == 999, ]
cordat$gdm <- factor(cordat$gdm)

names(cordat)[31] <- "birth_weight"
cordat$birth_weight <- scale(cordat$birth_weight)

# Filter to specific sequencing batches for quality control
cordat <- cordat[cordat$sequencing_batch %in% c("SK-3KKV", "SK-3KKW"), ]

# Prepare file paths for Gen3G inferential variance analysis
coldata.x <- cordat
names(coldata.x)[1] <- "names"
coldata.x$files <- paste0(DIR_COUNTS_assembly, "/", coldata.x$names, "/quant.sf")

coldata.y <- cordat
names(coldata.y)[1] <- "names"
coldata.y$files <- paste0(DIR_COUNTS_gencode, "/", coldata.y$names, "/quant.sf")

## Gen3G Assembly inferential variance
setTximetaBFC("ESPRESSO_corrected_SC_filtered")
makeLinkedTxome(indexDir = "ESPRESSO_corrected_SC_filtered",
                source = "ESPRESSO assembly of placenta lr-rnaseq, redo (filtered)",
                organism = "Homo sapiens",
                genome = "GRCh38",
                release = "p14",
                fasta = "ESPRESSO_corrected_SC_filtered/salmon_index/gentrome.fa",
                gtf = "ESPRESSO_corrected_SC_filtered.gtf",
                write = FALSE)

se.x <- tximeta(coldata.x[, c(32, 1:31)], type = "salmon")
se.x <- labelKeep(se.x)
se.x <- se.x[mcols(se.x)$keep, ]

infrv.assembly.Gen3G <- computeInfRV(se.x)
infrv.assembly.Gen3G <- mcols(infrv.assembly.Gen3G)$meanInfRV
rm(se.x)
gc()

## Gen3G GENCODE inferential variance
setTximetaBFC("GENCODE")
makeLinkedTxome(indexDir = "GENCODE/gencode.v45.salmon_index/gencode_v45",
                source = "GENCODEv45",
                organism = "Homo sapiens",
                genome = "GRCh38",
                release = "p14",
                fasta = "GENCODE/gencode.v45.salmon_index/gentrome.fa",
                gtf = "gencode_v45_annotation.gtf",
                write = FALSE)

se.y <- tximeta(coldata.y[, c(32, 1:31)], type = "salmon")
se.y <- labelKeep(se.y)
se.y <- se.y[mcols(se.y)$keep, ]

infrv.gencode.Gen3G <- computeInfRV(se.y)
infrv.gencode.Gen3G <- mcols(infrv.gencode.Gen3G)$meanInfRV
rm(se.y)
gc()

## Process and visualize Gen3G inferential uncertainty
infrv.Gen3G <- data.frame(transcript_id = union(names(infrv.assembly.Gen3G),
                                                names(infrv.gencode.Gen3G)))

infrv.assembly.Gen3G <- data.frame(infrv.assembly.Gen3G)
infrv.assembly.Gen3G$transcript_id <- row.names(infrv.assembly.Gen3G)
infrv.assembly.Gen3G <- infrv.assembly.Gen3G[, c(2, 1)]
names(infrv.assembly.Gen3G)[2] <- "lr-assembly"
infrv.Gen3G <- left_join(infrv.Gen3G, infrv.assembly.Gen3G)

infrv.gencode.Gen3G <- data.frame(infrv.gencode.Gen3G)
infrv.gencode.Gen3G$transcript_id <- row.names(infrv.gencode.Gen3G)
infrv.gencode.Gen3G <- infrv.gencode.Gen3G[, c(2, 1)]
names(infrv.gencode.Gen3G)[2] <- "GENCODEv45"
infrv.Gen3G <- left_join(infrv.Gen3G, infrv.gencode.Gen3G)

infrv.Gen3G <- infrv.Gen3G %>%
  pivot_longer(
    cols = c(`lr-assembly`, GENCODEv45),
    names_to = "annotation",
    values_to = "InfRV"
  )

infrv.Gen3G <- infrv.Gen3G[complete.cases(infrv.Gen3G), ]
names(infrv.Gen3G)[1] <- "isoform"
infrv.Gen3G <- left_join(infrv.Gen3G, classif.pcorr[, c("isoform", "structural_category")])
infrv.Gen3G$structural_category <- as.character(infrv.Gen3G$structural_category)
infrv.Gen3G[infrv.Gen3G$annotation == "GENCODEv45", "structural_category"] <- "GENCODE"
infrv.Gen3G$structural_category <- factor(infrv.Gen3G$structural_category,
                                          levels = c("GENCODE", "FSM", "ISM", "NIC", "NNC", "Other"))

## Figure S4D: Inferential uncertainty distribution (Gen3G)
figs4d <- ggplot(infrv.Gen3G, 
                 aes(x = log10(InfRV), fill = annotation, y = after_stat(density))) +
  geom_histogram(color = "#e9ecef", alpha = 0.4, position = 'identity', bins = 60) +
  scale_fill_manual(values = c("royalblue", "seagreen")) +
  theme_classic() +
  labs(fill = "") +
  ggtitle("") +
  geom_vline(xintercept = mean(log10(as.numeric(unlist(infrv.Gen3G[infrv.Gen3G$annotation == "lr-assembly", "InfRV"])))),
             color = "seagreen", linewidth = 1, linetype = "dashed") +
  geom_vline(xintercept = mean(log10(as.numeric(unlist(infrv.Gen3G[infrv.Gen3G$annotation == "GENCODEv45", "InfRV"])))),
             color = "royalblue", linewidth = 1, linetype = "dashed") +
  theme(panel.grid.major.y = element_line(),
        panel.grid.minor.y = element_line(),
        panel.grid.major.x = element_line(),
        panel.grid.minor.x = element_line(),
        legend.position = "top") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("mean InfRV") + 
  ylab("Density") +
  ggtitle("Isoform-level inferential uncertainty\nAll isoforms (Gen3G n=152)")

figs4d
ggsave("figs4d.png", figs4d, width = 6, height = 6, units = "in", dpi = 300)

## Figure S4E: Inferential uncertainty by structural category (Gen3G)
figs4e <- ggplot(infrv.Gen3G, aes(x = structural_category, fill = structural_category, y = log10(InfRV))) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.5, width = 0.1) +
  stat_halfeye(adjust = .5, width = .6, justification = -.2, .width = 0, point_colour = NA) +
  theme_minimal() + 
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.title.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_fill_manual(values = c("goldenrod", myPalette2)) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("mean InfRV") +
  ggtitle("Isoform-level inferential uncertainty\n(Gen3G n=152)")

figs4e
ggsave("figs4e.png", figs4e, width = 6, height = 6, units = "in", dpi = 300)

################################################################################
# Pairwise Sample Correlation Analysis
################################################################################

## Parallel correlation function for efficient computation
pcorr_parallel <- function(datExpr, num_cores = parallel::detectCores() - 1) {
  # Load required packages
  suppressPackageStartupMessages({
    library(parallel)
    library(foreach)
    library(doParallel)
  })
  
  # Register parallel backend
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Get column names and number of samples
  sample_names <- colnames(datExpr)
  num_samples <- ncol(datExpr)
  
  # Create all pairs of sample indices
  pairs <- expand.grid(i = 1:(num_samples-1), 
                       j = 2:num_samples)
  # Keep only pairs where j > i
  pairs <- pairs[pairs$j > pairs$i, ]
  
  # Process in parallel
  corr_results <- foreach(idx = 1:nrow(pairs), .combine = rbind) %dopar% {
    i <- pairs$i[idx]
    j <- pairs$j[idx]
    sample_A <- sample_names[i]
    sample_B <- sample_names[j]
    corr_test <- suppressWarnings(cor.test(datExpr[, i], datExpr[, j], method = "pearson"))
    c(sample_A, sample_B, corr_test$estimate)
  }
  
  # Stop cluster
  stopCluster(cl)
  
  # Create dataframe and format results
  exprCorr <- data.frame(corr_results)
  colnames(exprCorr) <- c("Sample_A", "Sample_B", "R")
  exprCorr$R <- as.numeric(exprCorr$R)
  
  return(exprCorr)
}

## Function to calculate correlations by structural category
pcorr.bystruc <- function(datExpr, tlens, classif, category) {
  # Convert to TPM-like values
  datExpr <- datExpr / tlens
  datExpr <- datExpr[rownames(datExpr) %in% classif[classif$structural_category == category, "transcript_id"], ]
  datExpr <- log(t(t(datExpr) * 1e6 / colSums(datExpr)) + 1)
  
  # Calculate pairwise correlations
  exprCorr <- pcorr_parallel(datExpr)
  exprCorr$group <- category
  return(exprCorr)
}

# Adjust classification naming for correlation analysis
names(classif.pcorr)[1] <- "transcript_id"

################################################################################
# GUSTO Pairwise Sample Correlations
################################################################################

## Calculate correlations by structural category
counts.GUSTO.assembly <- counts(estimateSizeFactors(dds.assembly.GUSTO), normalized = TRUE)
GUSTO.exprCorr.assembly.FSM <- pcorr.bystruc(counts.GUSTO.assembly, txi.assembly.GUSTO$length, classif.pcorr, "FSM")
GUSTO.exprCorr.assembly.ISM <- pcorr.bystruc(counts.GUSTO.assembly, txi.assembly.GUSTO$length, classif.pcorr, "ISM")
GUSTO.exprCorr.assembly.NIC <- pcorr.bystruc(counts.GUSTO.assembly, txi.assembly.GUSTO$length, classif.pcorr, "NIC")
GUSTO.exprCorr.assembly.NNC <- pcorr.bystruc(counts.GUSTO.assembly, txi.assembly.GUSTO$length, classif.pcorr, "NNC")
GUSTO.exprCorr.assembly.Other <- pcorr.bystruc(counts.GUSTO.assembly, txi.assembly.GUSTO$length, classif.pcorr, "Other")

## GENCODE correlations for comparison
counts.GUSTO.gencode <- counts(estimateSizeFactors(dds.gencode.GUSTO), normalized = TRUE)
GUSTO.datExpr.TPM.gencode <- counts.GUSTO.gencode / txi.gencode.GUSTO$length
GUSTO.datExpr.TPM.gencode <- log(t(t(GUSTO.datExpr.TPM.gencode) * 1e6 / colSums(GUSTO.datExpr.TPM.gencode)) + 1)
GUSTO.exprCorr.gencode <- pcorr_parallel(GUSTO.datExpr.TPM.gencode)
GUSTO.exprCorr.gencode$group <- "GENCODE"

# Combine all correlation results
GUSTO.exprCorr.merged <- rbind(GUSTO.exprCorr.assembly.FSM, GUSTO.exprCorr.assembly.ISM,
                               GUSTO.exprCorr.assembly.NIC, GUSTO.exprCorr.assembly.NNC,
                               GUSTO.exprCorr.assembly.Other, GUSTO.exprCorr.gencode)
GUSTO.exprCorr.merged$group <- factor(GUSTO.exprCorr.merged$group,
                                      levels = c("GENCODE", "FSM", "ISM", "NIC", "NNC", "Other"))

## Figure S4F: GUSTO pairwise sample correlations
figs4f <- ggplot(GUSTO.exprCorr.merged, aes(x = group, fill = group, y = R)) +
  geom_boxplot(alpha = 0.5, outliers = FALSE) +
  theme_minimal() + 
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.title.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_fill_manual(values = c("goldenrod", myPalette2)) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("R") + 
  ylim(0.5, 1) +
  ggtitle("Pairwise Sample Correlations\n(GUSTO n=200)")

figs4f
ggsave("figs4f.png", figs4f, width = 5, height = 5, units = "in", dpi = 300)

################################################################################
# Gen3G Pairwise Sample Correlations
################################################################################

## Calculate correlations by structural category for Gen3G
counts.Gen3G.assembly <- counts(estimateSizeFactors(dds.assembly.Gen3G), normalized = TRUE)
Gen3G.exprCorr.assembly.FSM <- pcorr.bystruc(counts.Gen3G.assembly, txi.assembly.Gen3G$length, classif.pcorr, "FSM")
Gen3G.exprCorr.assembly.ISM <- pcorr.bystruc(counts.Gen3G.assembly, txi.assembly.Gen3G$length, classif.pcorr, "ISM")
Gen3G.exprCorr.assembly.NIC <- pcorr.bystruc(counts.Gen3G.assembly, txi.assembly.Gen3G$length, classif.pcorr, "NIC")
Gen3G.exprCorr.assembly.NNC <- pcorr.bystruc(counts.Gen3G.assembly, txi.assembly.Gen3G$length, classif.pcorr, "NNC")
Gen3G.exprCorr.assembly.Other <- pcorr.bystruc(counts.Gen3G.assembly, txi.assembly.Gen3G$length, classif.pcorr, "Other")

## Gen3G GENCODE correlations
counts.Gen3G.gencode <- counts(estimateSizeFactors(dds.gencode.Gen3G), normalized = TRUE)
Gen3G.datExpr.TPM.gencode <- counts.Gen3G.gencode / txi.gencode.Gen3G$length
Gen3G.datExpr.TPM.gencode <- log(t(t(Gen3G.datExpr.TPM.gencode) * 1e6 / colSums(Gen3G.datExpr.TPM.gencode)) + 1)
Gen3G.exprCorr.gencode <- pcorr_parallel(Gen3G.datExpr.TPM.gencode)
Gen3G.exprCorr.gencode$group <- "GENCODE"

# Combine Gen3G correlation results
Gen3G.exprCorr.merged <- rbind(Gen3G.exprCorr.assembly.FSM, Gen3G.exprCorr.assembly.ISM,
                               Gen3G.exprCorr.assembly.NIC, Gen3G.exprCorr.assembly.NNC,
                               Gen3G.exprCorr.assembly.Other, Gen3G.exprCorr.gencode)
Gen3G.exprCorr.merged$group <- factor(Gen3G.exprCorr.merged$group,
                                      levels = c("GENCODE", "FSM", "ISM", "NIC", "NNC", "Other"))

## Figure S4G: Gen3G pairwise sample correlations
figs4g <- ggplot(Gen3G.exprCorr.merged, aes(x = group, fill = group, y = R)) +
  geom_boxplot(alpha = 0.5, outliers = FALSE) +
  theme_minimal() + 
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.title.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_fill_manual(values = c("goldenrod", myPalette2)) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("R") + 
  ylim(0.5, 1) +
  ggtitle("Pairwise Sample Correlations\n(Gen3G n=152)")

figs4g
ggsave("figs4g.png", figs4g, width = 5, height = 5, units = "in", dpi = 300)

################################################################################
# Transcript Feature Analysis
################################################################################

## Classify transcripts by novelty for feature analysis
novelGtx <- classif.filtered[grep("novel", classif.filtered$associated_gene), "isoform"]
fsm <- classif.filtered[classif.filtered$structural_category %in% c("FSM", "ISM"), "isoform"]
noveltx <- setdiff(classif.filtered[!classif.filtered$structural_category %in% c("FSM", "ISM"), "isoform"],
                   novelGtx)

################################################################################
# Figure S5A: Transcript abundance
################################################################################

fl.dat <- data.frame(FL=c(classif.filtered[classif.filtered$isoform%in%fsm,"FL"],
                          classif.filtered[classif.filtered$isoform%in%noveltx,"FL"],
                          classif.filtered[classif.filtered$isoform%in%c(fsm,noveltx),"FL"],
                          classif.filtered[classif.filtered$isoform%in%novelGtx,"FL"]),
                     set=c(rep.int("Known",length(classif.filtered[classif.filtered$isoform%in%fsm,"FL"])),
                           rep.int("Novel",length(classif.filtered[classif.filtered$isoform%in%noveltx,"FL"])),
                           rep.int("Known",length(classif.filtered[classif.filtered$isoform%in%c(fsm,noveltx),"FL"])),
                           rep.int("Novel",length(classif.filtered[classif.filtered$isoform%in%novelGtx,"FL"]))),
                     comparison=c(rep.int("Within known genes",length(classif.filtered[classif.filtered$isoform%in%fsm,"FL"])+length(classif.filtered[classif.filtered$isoform%in%noveltx,"FL"])),
                                  rep.int("Known vs novel genes",length(classif.filtered[classif.filtered$isoform%in%c(fsm,noveltx),"FL"])+length(classif.filtered[classif.filtered$isoform%in%novelGtx,"FL"]))))

figs5a <- ggplot(fl.dat, aes(x = comparison, y = log10(FL), color = set)) +
  geom_boxplot() + 
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Full-length reads (log10)") +
  scale_color_manual(values = c("black", "grey"),
                     name = "Transcripts")

figs5a
ggsave("figs5a.png", figs5a, width = 4, height = 4, units = "in", dpi = 300)

# Statistical testing
wilcox.test(classif.filtered[classif.filtered$isoform %in% fsm, "FL"],
            classif.filtered[classif.filtered$isoform %in% noveltx, "FL"],
            alternative = "greater")

wilcox.test(classif.filtered[classif.filtered$isoform %in% c(fsm, noveltx), "FL"],
            classif.filtered[classif.filtered$isoform %in% novelGtx, "FL"],
            alternative = "greater")

################################################################################
# Figure S5B: Transcript Length by Category
################################################################################

len.dat <- data.frame(length = c(classif.filtered[classif.filtered$isoform %in% fsm, "length"],
                                 classif.filtered[classif.filtered$isoform %in% noveltx, "length"],
                                 classif.filtered[classif.filtered$isoform %in% c(fsm, noveltx), "length"],
                                 classif.filtered[classif.filtered$isoform %in% novelGtx, "length"]),
                      set = c(rep.int("Known", length(classif.filtered[classif.filtered$isoform %in% fsm, "length"])),
                              rep.int("Novel", length(classif.filtered[classif.filtered$isoform %in% noveltx, "length"])),
                              rep.int("Known", length(classif.filtered[classif.filtered$isoform %in% c(fsm, noveltx), "length"])),
                              rep.int("Novel", length(classif.filtered[classif.filtered$isoform %in% novelGtx, "length"]))),
                      comparison = c(rep.int("Within known genes", length(classif.filtered[classif.filtered$isoform %in% fsm, "length"]) + length(classif.filtered[classif.filtered$isoform %in% noveltx, "length"])),
                                     rep.int("Known vs novel genes", length(classif.filtered[classif.filtered$isoform %in% c(fsm, noveltx), "length"]) + length(classif.filtered[classif.filtered$isoform %in% novelGtx, "length"]))))

figs5b <- ggplot(len.dat, aes(x = comparison, y = log10(length), color = set)) +
  geom_boxplot() + 
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Transcript length (log10 bp)") +
  scale_color_manual(values = c("black", "grey"),
                     name = "Transcripts")

figs5b
ggsave("figs5b.png", figs5b, width = 4, height = 4, units = "in", dpi = 300)

# Statistical testing
wilcox.test(classif.filtered[classif.filtered$isoform %in% fsm, "length"],
            classif.filtered[classif.filtered$isoform %in% noveltx, "length"],
            alternative = "greater")

wilcox.test(classif.filtered[classif.filtered$isoform %in% c(fsm, noveltx), "length"],
            classif.filtered[classif.filtered$isoform %in% novelGtx, "length"],
            alternative = "greater")

################################################################################
# Figure S5C: Exon Count by Category
################################################################################

exons.dat <- data.frame(length = c(classif.filtered[classif.filtered$isoform %in% fsm, "exons"],
                                   classif.filtered[classif.filtered$isoform %in% noveltx, "exons"],
                                   classif.filtered[classif.filtered$isoform %in% c(fsm, noveltx), "exons"],
                                   classif.filtered[classif.filtered$isoform %in% novelGtx, "exons"]),
                        set = c(rep.int("Known", length(classif.filtered[classif.filtered$isoform %in% fsm, "exons"])),
                                rep.int("Novel", length(classif.filtered[classif.filtered$isoform %in% noveltx, "exons"])),
                                rep.int("Known", length(classif.filtered[classif.filtered$isoform %in% c(fsm, noveltx), "exons"])),
                                rep.int("Novel", length(classif.filtered[classif.filtered$isoform %in% novelGtx, "exons"]))),
                        comparison = c(rep.int("Within known genes", length(classif.filtered[classif.filtered$isoform %in% fsm, "exons"]) + length(classif.filtered[classif.filtered$isoform %in% noveltx, "exons"])),
                                       rep.int("Known vs novel genes", length(classif.filtered[classif.filtered$isoform %in% c(fsm, noveltx), "exons"]) + length(classif.filtered[classif.filtered$isoform %in% novelGtx, "exons"]))))

figs5c <- ggplot(exons.dat, aes(x = comparison, y = log10(length), color = set)) +
  geom_boxplot() + 
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Exons per transcript") +
  scale_color_manual(values = c("black", "grey"),
                     name = "Transcripts")

figs5c
ggsave("figs5c.png", figs5c, width = 4, height = 4, units = "in", dpi = 300)

# Statistical testing
wilcox.test(classif.filtered[classif.filtered$isoform %in% fsm, "exons"],
            classif.filtered[classif.filtered$isoform %in% noveltx, "exons"],
            alternative = "greater")

wilcox.test(classif.filtered[classif.filtered$isoform %in% c(fsm, noveltx), "exons"],
            classif.filtered[classif.filtered$isoform %in% novelGtx, "exons"],
            alternative = "greater")

################################################################################
# Coding Potential Analysis
################################################################################

## Import coding potential analysis results
CPC2 <- read.table("/Users/stbresnahan/Desktop/Manuscripts/Placenta/CPC2.txt", header = TRUE)
names(CPC2)[1] <- "isoform"
CPC2$isoform <- sub("\\.[^.]*$", "", CPC2$isoform)
CPC2 <- left_join(CPC2, classif.filtered[, c("isoform", "structural_category")])
CPC2 <- CPC2[!is.na(CPC2$structural_category), ]

## Import ORF prediction results
TD_ORFs <- read.table("/Users/stbresnahan/Desktop/Manuscripts/Placenta/longest_orfs.gff3", header = FALSE)
TD_ORFs <- TD_ORFs[TD_ORFs$V3 == "CDS", c(1, 5)]
TD_ORFs <- TD_ORFs[!duplicated(TD_ORFs), ]
names(TD_ORFs) <- c("isoform", "ORF_length")

TD_ORFs <- TD_ORFs %>%
  group_by(isoform) %>%
  slice_max(order_by = ORF_length, n = 1, with_ties = FALSE) %>%
  ungroup()

TD_ORFs <- left_join(TD_ORFs, classif.filtered[, c("isoform", "structural_category")])
TD_ORFs <- TD_ORFs[!is.na(TD_ORFs$structural_category), ]

################################################################################
# Figure S6A: Percentage of Isoforms with ORFs by Category
################################################################################

classif.filtered$ORF <- TRUE
classif.filtered[is.na(classif.filtered$ORF_seq), "ORF"] <- FALSE
ORF.td <- data.frame(table(classif.filtered[, c("structural_category", "ORF")]))

temp <- ORF.td %>%
  mutate(ORF = as.character(ORF)) %>%
  pivot_wider(names_from = ORF, values_from = Freq, values_fill = 0)

ORF_percent <- data.frame(temp)
names(ORF_percent) <- c("structural_category", "FALSE", "TRUE")
ORF_percent$percent_T <- ORF_percent$`TRUE` / sum(ORF_percent$`FALSE` + ORF_percent$`TRUE`)

figs6a <- ggplot(ORF_percent, aes(x = structural_category, y = percent_T,
                                  fill = structural_category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  scale_fill_manual(values = myPalette, guide = 'none') +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.title.x = element_blank()) +
  labs(y = "% of isoforms with ORFs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

figs6a
ggsave("figs6a.png", figs6a, width = 4, height = 4, units = "in", dpi = 300)

################################################################################
# Figure S6B: ORF Length by Structural Category
################################################################################

figs6b <- ggplot(TD_ORFs, aes(x = structural_category, y = log10(ORF_length + 1), fill = structural_category)) +
  geom_boxplot(outliers = FALSE) + 
  scale_fill_manual(values = myPalette) + 
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  ylab("ORF length (log10)")

figs6b
ggsave("figs6b.png", figs6b, width = 4, height = 4, units = "in", dpi = 300)

################################################################################
# Figure S6C: Coding Probability by Structural Category
################################################################################

figs6c <- ggplot(CPC2, aes(x = structural_category, y = coding_probability, fill = structural_category)) +
  geom_boxplot(outliers = FALSE) + 
  scale_fill_manual(values = myPalette) + 
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  ylab("Coding Probability")

figs6c
ggsave("figs6c.png", figs6c, width = 4, height = 4, units = "in", dpi = 300)

################################################################################
# CDS Support Analysis
################################################################################

## BLASTP analysis with placenta-specific peptides
BLASTP <- read.table("/Users/stbresnahan/Desktop/Manuscripts/Placenta/ESPRESSO_corrected_SC_filtered/ESPRESSO_corrected_SC_filtered_blastp.outfmt6", header = FALSE)
BLASTP <- BLASTP[, c(1, 2, 11, 12)]
names(BLASTP) <- c("isoform", "protein_accession", "Evalue", "score")
BLASTP <- left_join(BLASTP, classif.filtered[, c("isoform", "structural_category")])
BLASTP <- BLASTP[!is.na(BLASTP$structural_category), ]

CDS_support <- intersect(BLASTP[BLASTP$Evalue < 0.01 & BLASTP$score > 10, "isoform"],
                         CPC2[CPC2$coding_probability > 0.5, "isoform"])

classif.filtered$CDS_support <- FALSE
classif.filtered[classif.filtered$isoform %in% CDS_support, "CDS_support"] <- TRUE

classif.filtered$TD_ORF <- FALSE
classif.filtered[classif.filtered$isoform %in% TD_ORFs$isoform, "TD_ORF"] <- TRUE

CStab <- data.frame(table(classif.filtered[, c("structural_category", "TD_ORF", "CDS_support")]))
CStab <- CStab %>%
  group_by(structural_category, TD_ORF, CDS_support) %>%
  pivot_wider(names_from = CDS_support, values_from = Freq, values_fill = 0)
CStab <- data.frame(CStab)
names(CStab) <- c("structural_category", "ORF", "False", "True")
CStab$percentT <- CStab$True / (CStab$False + CStab$True)
CStab <- CStab[10:18, ]

################################################################################
# Figure 2C: CDS Support from Placenta MS/MS
################################################################################

fig2c <- ggplot(CStab, aes(x = structural_category, y = percentT, fill = structural_category)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = myPalette) +
  theme_minimal() +
  theme(axis.title.x = element_blank(), legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Coding isoforms with\nCDS support",
       y = "% of coding isoforms") +
  scale_y_continuous(labels = scales::percent)

fig2c
ggsave("fig2c.png", fig2c, dpi = 300, width = 3, height = 3, units = "in")

## BLASTP analysis with full UniProt database
BLASTP.full <- read.table("/Users/stbresnahan/Desktop/Manuscripts/Placenta/blastp.outfmt6", header = FALSE)
BLASTP.full <- BLASTP.full[, c(1, 2, 11, 12)]
names(BLASTP.full) <- c("isoform", "protein_accession", "Evalue", "score")
BLASTP.full$protein_accession <- sub("^[^|]*\\|([^|]*)\\|.*$", "\\1", BLASTP.full$protein_accession)
BLASTP.full$isoform <- sub("\\.[^.]*$", "", BLASTP.full$isoform)
BLASTP.full <- left_join(BLASTP.full, classif.filtered[, c("isoform", "structural_category")])
BLASTP.full <- BLASTP.full[!is.na(BLASTP.full$structural_category), ]

CDS_support_full <- intersect(BLASTP.full[BLASTP.full$Evalue < 0.01 & BLASTP.full$score > 10, "isoform"],
                              CPC2[CPC2$coding_probability > 0.5, "isoform"])

classif.filtered$CDS_support_full <- FALSE
classif.filtered[classif.filtered$isoform %in% CDS_support_full, "CDS_support_full"] <- TRUE

CStab <- data.frame(table(classif.filtered[, c("structural_category", "TD_ORF", "CDS_support_full")]))
CStab <- CStab %>%
  group_by(structural_category, TD_ORF, CDS_support_full) %>%
  pivot_wider(names_from = CDS_support_full, values_from = Freq, values_fill = 0)
CStab <- data.frame(CStab)
names(CStab) <- c("structural_category", "ORF", "False", "True")
CStab$percentT <- CStab$True / (CStab$False + CStab$True)
CStab <- CStab[10:18, ]

################################################################################
# Figure S6D: CDS Support from Full UniProt
################################################################################

figs6d <- ggplot(CStab, aes(x = structural_category, y = percentT, fill = structural_category)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = myPalette) +
  theme_minimal() +
  theme(axis.title.x = element_blank(), legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Coding isoforms with\nCDS support",
       y = "% of coding isoforms") +
  scale_y_continuous(labels = scales::percent)

figs6d
ggsave("figs6d.png", figs6d, dpi = 300, width = 3, height = 3, units = "in")

################################################################################
# Gene Enrichment Analysis
################################################################################

## GO enrichment of top 50 genes with most isoforms
citab <- data.frame(table(classif.filtered$associated_gene))
names(citab)[1] <- "gene"
citab <- citab[rev(order(citab$Freq)), ]

tx2g.assembly <- classif.filtered[, c("isoform", "associated_gene")]
tx2g.assembly$isoform <- sub("\\.[^.]*$", "", tx2g.assembly$isoform)

## Function to convert Ensembl IDs to gene symbols
ensembl_to_symbol <- function(ensembl_ids) {
  ensembl_ids <- sapply(strsplit(ensembl_ids, '[.]'), function(x) x[1])
  
  gene_symbols <- mapIds(
    org.Hs.eg.db,
    keys = ensembl_ids,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  # Remove NAs
  valid_idx <- !is.na(gene_symbols)
  genes <- gene_symbols[valid_idx]
  ensembl <- names(genes)
  
  mapping <- data.frame(
    ensembl_id = ensembl,
    gene_symbol = genes,
    stringsAsFactors = FALSE
  )
  
  return(list(symbols = as.character(genes), mapping = mapping))
}

symbols <- ensembl_to_symbol(unique(as.character(citab[1:50, "gene"])))

setEnrichrSite("Enrichr")
dbs <- c(
  "GO_Biological_Process_2023",
  "GO_Molecular_Function_2023", 
  "GO_Cellular_Component_2023",
  "KEGG_2021_Human",
  "WikiPathway_2021_Human",
  "Reactome_2022",
  "MSigDB_Hallmark_2020"
)

enrichr_results <- enrichr(
  genes = symbols$symbols,
  databases = dbs,
  background = ensembl_to_symbol(unique(as.character(tx2g.assembly$associated_gene)))$symbols
)

print(enrichr_results[["GO_Molecular_Function_2023"]][["Term"]][1])
## > "Hormone Receptor Binding (GO:0051427)"

################################################################################
# Transcript-Gene Correlation Analysis
################################################################################

## Correlations of isoforms per gene with various features
citab <- data.frame(table(classif.filtered$associated_gene))
names(citab) <- c("associated_gene", "isoforms")
max_lengths <- aggregate(length ~ associated_gene, data = classif.filtered, FUN = max)
sum_exons <- aggregate(exons ~ associated_gene, data = classif.filtered, FUN = sum)
sum_FL <- aggregate(FL ~ associated_gene, data = classif.filtered, FUN = sum)
citab <- left_join(citab, max_lengths)
citab <- left_join(citab, sum_exons)
citab <- left_join(citab, sum_FL)

classif.filtered$RPK <- classif.filtered$FL / (classif.filtered$length / 1000)
classif.filtered$TPM <- (classif.filtered$RPK / sum(classif.filtered$RPK)) * 1e6
sum_TPM <- aggregate(TPM ~ associated_gene, data = classif.filtered, FUN = sum)
citab <- left_join(citab, sum_TPM)
citab$log10TPM <- log10(citab$TPM)

################################################################################
# Figure S5D: Transcripts per Gene vs Exons (All Transcripts)
################################################################################

t2 <- cor.test(x = citab$isoforms, y = citab$exons)
r_val <- round(t2$estimate, 3)
p_val <- signif(t2$p.value, 3)
label_text <- paste0("r = ", r_val, "\np = ", p_val)

exoncor.df <- data.frame(isoforms = citab$isoforms,
                         exons = citab$exons)

figs5d <- ggplot(exoncor.df, aes(x = isoforms, y = exons)) +
  geom_point() +  
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal() +
  xlab("Transcripts per gene") +
  ylab("Exons per transcript") +
  annotate("text", x = 30, y = 500, label = label_text, hjust = 0) +
  ggtitle("All transcripts")

figs5d
ggsave("figs5d.png", figs5d, width = 4, height = 4, units = "in", dpi = 300)

################################################################################
# Figure S5E: Transcripts per Gene vs Exons (Highly-Expressed Transcripts)
################################################################################

t4 <- cor.test(x = citab[citab$log10TPM > 2.5, ]$isoforms, y = citab[citab$log10TPM > 2.5, ]$exons)
r_val <- round(t4$estimate, 3)
p_val <- signif(t4$p.value, 3)
label_text <- paste0("r = ", r_val, "\np = ", p_val)

exoncor.df <- data.frame(isoforms = citab[citab$log10TPM > 2.5, ]$isoforms,
                         exons = citab[citab$log10TPM > 2.5, ]$exons)

figs5e <- ggplot(exoncor.df, aes(x = isoforms, y = exons)) +
  geom_point() +  
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal() +
  xlab("Transcripts per gene") +
  ylab("Exons per transcript") +
  annotate("text", x = 30, y = 500, label = label_text, hjust = 0) +
  ggtitle("log10 TPM > 2.5")

figs5e
ggsave("figs5e.png", figs5e, width = 4, height = 4, units = "in", dpi = 300)

################################################################################
# End of Analysis
################################################################################