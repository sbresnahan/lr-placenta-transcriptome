################################################################################
# Short-Read Transcriptomics Analysis: Assembly vs Reference Comparison
# Author: Sean T. Bresnahan
# Description: Comprehensive comparison of short-read RNA-seq quantification using
#              the filtered placental transcriptome assembly versus GENCODE v45
#              reference. Analyzes expression patterns, statistical properties,
#              quantification uncertainty, transcript features, and coding potential
#              across two independent cohorts (GUSTO N=200, Gen3G N=152).
#
# Input Datasets:
#   - Assembly/ESPRESSO_corrected_SC_filtered_classification.txt
#       SQANTI3 structural classification for filtered assembly
#   - Assembly/ESPRESSO_corrected_classification.txt
#       SQANTI3 structural classification for unfiltered assembly
#   - Assembly/ESPRESSO_corrected_SC_filtered.gtf
#       GTF annotation for filtered assembly
#   - Assembly/ESPRESSO_corrected_corrected.gtf
#       GTF annotation for unfiltered assembly (for filtering comparison)
#   - Assembly/lr_assembly_GENCODE_v45_combined.gtf
#       Combined GENCODE v45 + filtered assembly GTF (GENCODE+)
#   - /rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v45/gencode.v45.annotation.gtf
#       GENCODE v45 reference annotation
#   - /rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/metadata_GUSTO.csv
#       GUSTO cohort metadata (N=200 samples) derived from GUSTO data vault documentation
#   - /rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/infrv_GUSTO/COUNTS_assembly_SC/*/quant.sf
#       Salmon quantification files for GUSTO cohort using lr-assembly
#   - /rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/infrv_GUSTO/COUNTS_gencode/*/quant.sf
#       Salmon quantification files for GUSTO cohort using GENCODE v45
#   - /rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/infrv_GUSTO/COUNTS_combined/*/quant.sf
#       Salmon quantification files for GUSTO cohort using GENCODE+
#   - /rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/metadata_Gen3G.csv
#       Gen3G cohort metadata files from dbGaP (N=152 samples) derived from dbGaP phenotype accessions
#   - /rsrch5/home/epi/bhattacharya_lab/data/mapqtl/Gen3G/metadata/SraRunTable.csv
#       Gen3G SRA run table with sample mappings
#   - /rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/infrv_Gen3G/COUNTS_assembly_SC/*/quant.sf
#       Salmon quantification files for Gen3G cohort using lr-assembly
#   - /rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/infrv_Gen3G/COUNTS_gencode/*/quant.sf
#       Salmon quantification files for Gen3G cohort using GENCODE v45
#   - /rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/infrv_Gen3G/COUNTS_combined/*/quant.sf
#       Salmon quantification files for Gen3G cohort using GENCODE+
#   - /rsrch5/home/epi/stbresnahan/scratch/natcomms_revision/COUNTS_GUSTO_Adipose_Subcutaneous/*/quant.sf
#       Salmon quantification for GUSTO samples using GTEx Adipose assembly
#   - /rsrch5/home/epi/stbresnahan/scratch/natcomms_revision/COUNTS_GUSTO_Cells_Cultured_fibroblasts/*/quant.sf
#       Salmon quantification for GUSTO samples using GTEx Fibroblasts assembly
#
# Output Files:
#   Data Files:
#     - Supplementary/tx2g_ESPRESSO_assembly_SC_filtered.RData
#         Transcript-to-gene mapping objects (tx2g.assembly, tx2g.gencode, tx2g.combined)
#     - Supplementary/metadata_GUSTO.RData
#         GUSTO cohort metadata objects (metadata.GUSTO, cordat.GUSTO)
#     - Supplementary/metadata_Gen3G.RData
#         Gen3G cohort metadata objects (metadata.Gen3G, cordat.Gen3G.Gen3G)
#
#   Figures - Sample Quality Control:
#     - Figures/tx_complex_by_sex.png
#         Transcriptional complexity (isoforms per gene) by fetal sex
#     - Figures/figs6a.pdf
#         PCA plots for GUSTO and Gen3G cohorts (quality control)
#
#   Figures - Expression Comparisons:
#     - Figures/fig3b.pdf
#         Total isoform expression across samples (GUSTO and Gen3G)
#     - Figures/fig3c.pdf
#         Mean and variance by transcript categories (GUSTO, discovery vs filtration)
#     - Figures/figs6b.pdf
#         Mean and variance by transcript categories (Gen3G, discovery vs filtration)
#
#   Figures - Inferential Uncertainty:
#     - Figures/fig3d.pdf
#         InfRV distribution by annotation (GUSTO: GENCODEv45, GENCODE+, lr-assembly)
#     - Figures/figs6c.pdf
#         InfRV distribution by annotation (Gen3G: GENCODEv45, GENCODE+, lr-assembly)
#     - Figures/figs6d.pdf
#         InfRV by structural category (GUSTO, lr-assembly only)
#     - Figures/figs6e.pdf
#         InfRV by structural category (Gen3G, lr-assembly only)
#     - Figures/figs6h.pdf
#         InfRV vs exon sharing scatter plot (GUSTO)
#     - Figures/figs6i.pdf
#         InfRV across placenta and GTEx tissues (GUSTO samples)
#
#   Figures - Sample Correlations:
#     - Figures/figs6f.pdf
#         Pairwise sample correlations by structural category (GUSTO)
#     - Figures/figs6g.pdf
#         Pairwise sample correlations by structural category (Gen3G)
#
################################################################################

################################################################################
# Environment Setup
################################################################################
# .libPaths(c( "/home/stbresnahan/R/ubuntu/4.3.1" , .libPaths()))

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
library(patchwork)     # Plot arrangement
library(doSNOW)
library(GenomicRanges)
library(rtracklayer)
library(broom)

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
classif <- read.table("Assembly/ESPRESSO_corrected_SC_filtered_classification.txt", header = TRUE)

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
classif.original <- read.table("Assembly/ESPRESSO_corrected_classification.txt", header = TRUE)
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
# Transcript-to-Gene Mapping for Filtered Assembly
################################################################################

# Import filtered assembly GTF and create transcript-to-gene mapping
tx2g.assembly <- data.frame(import("Assembly/ESPRESSO_corrected_SC_filtered.gtf"))
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
# GUSTO Cohort Analysis
################################################################################

# Prepare metadata
cordat.GUSTO <- read.csv("/rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/metadata_GUSTO.csv") # From GUSTO data vault
cordat.GUSTO <- metadata.GUSTO
cordat.GUSTO$GROUP <- factor(as.numeric(factor(cordat.GUSTO$GROUP)))
cordat.GUSTO$sex <- factor(as.numeric(factor(cordat.GUSTO$sex)))
cordat.GUSTO$gdm <- factor(cordat.GUSTO$gdm)
cordat.GUSTO$mother_ethnicity <- factor(cordat.GUSTO$mother_ethnicity)
cordat.GUSTO$htn_clean <- factor(as.numeric(factor(cordat.GUSTO$htn_clean)))
cordat.GUSTO$GA <- as.numeric(scale(cordat.GUSTO$GA))
cordat.GUSTO$c_weight_birth_kg <- as.numeric(scale(cordat.GUSTO$c_weight_birth_kg*1000))
names(cordat.GUSTO)[which(names(cordat.GUSTO)=="c_weight_birth_kg")] <- "birth_weight"
cordat <- cordat.GUSTO

## Prepare assembly gene expression data
dirs.salmon.assembly <- list.dirs("/rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/infrv_GUSTO/COUNTS_assembly_SC", recursive = FALSE)
files.list.assembly <- list()
for (i in 1:length(dirs.salmon.assembly)) {
  files.list.assembly[i] <- paste(unlist(dirs.salmon.assembly[i]), "quant.sf", sep = "/")
}
dirnames.assembly <- sapply(strsplit(unlist(dirs.salmon.assembly), "/"), "[",
                            length(unlist(strsplit(unlist(dirs.salmon.assembly), "/")[1])))
names(files.list.assembly) <- dirnames.assembly

# Import gene-level quantification
txi.assembly.g <- tximport(unlist(files.list.assembly), type = "salmon", tx2gene = tx2g.assembly,
                           geneIdCol = "gene_id", txIdCol = "transcript_id", dropInfReps = TRUE, txOut = FALSE)

### Quality control and outlier removal
dds.assembly <- DESeqDataSetFromTximport(txi.assembly.g, metadata, ~1)
vst.assembly <- vst(dds.assembly, blind = TRUE)
pca.assembly.check <- plotPCA(vst.assembly, intgroup = "sex")
pca.assembly.check

# Identify outliers from PCA
pca.assembly.dat <- pca.assembly.check[["data"]]
outlier_samples <- pca.assembly.dat[pca.assembly.dat$PC2 > 60, "name"] # "J1102" "J1147"

# Define samples to filter (outliers + missing GDM data)
filter.assembly <- c("J1102", "J1147", "J1106", "J1147",
                     metadata[is.na(metadata$gdm), "SampleID"])
txi.assembly.filter <- filter_tximport(txi.assembly.g, filter.assembly)

# Update metadata after filtering
cordat <- cordat[!cordat$SampleID %in% filter.assembly, ]

# Save metadata for later use
metadata.GUSTO <- metadata
cordat.GUSTO <- cordat
save(list=c("cordat.GUSTO"),file="Supplementary/metadata_GUSTO.RData")

################################################################################
# GUSTO Assembly Analysis
################################################################################

## Import Salmon quantification results and apply QU (quasi-UMI) correction
dirs.salmon.assembly <- list.dirs("/rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/infrv_GUSTO/COUNTS_assembly_SC", recursive = FALSE)

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

# PCA of GUSTO samples for outlier detection
ggplot(pca.df, aes(x = PC1, y = PC2, color = gdm)) +
  geom_point() +
  xlab("PC1: 25.72% variance") +
  ylab("PC2: 19.26% variance") +
  ggtitle("GUSTO Sample PCA")


##### Transcriptional complexity by fetal sex
counts <- data.frame(counts.GUSTO.assembly)
counts.M <- counts[,names(counts)%in%cordat.GUSTO[cordat.GUSTO$sex==1,"SampleID"]]
counts.F <- counts[,names(counts)%in%cordat.GUSTO[cordat.GUSTO$sex==2,"SampleID"]]

counts.M <- counts.M[rowSums(counts.M)>10,]
counts.F <- counts.F[rowSums(counts.F)>10,]

tab.M <- data.frame(table(tx2g.assembly[tx2g.assembly$transcript_id%in%row.names(counts.M),"gene_id"]))
tab.F <- data.frame(table(tx2g.assembly[tx2g.assembly$transcript_id%in%row.names(counts.F),"gene_id"]))

tab.res <- data.frame(freq=c(tab.M$Freq,tab.F$Freq),
                      group=c(rep.int("M",length(tab.M$Freq)),
                              rep.int("F",length(tab.F$Freq))))

t_test_result <- t.test(log(freq) ~ group, data = tab.res)
label_text <- paste0("Difference: ", round(diff(t_test_result$estimate), 2), 
                     "\np = ", round(t_test_result$p.value, 3))

g1 <- ggplot(tab.res,aes(x=group,y=log(freq))) +
  geom_boxplot() +
  annotate("text", x = 1.5, y = max(log(tab.res$freq), na.rm = TRUE), 
           label = label_text, vjust = 1) +
  labs(x="Fetal sex",y="Isoforms per gene (log)",
       title="Transcriptional complexity by fetal sex",
       subtitle="Data from GUSTO (n=200)") +
  theme_minimal()
g1

ggsave("Figures/tx_complex_by_sex.png",g1,width=6.57,height=3.85,units="in",dpi=300)

################################################################################
# GUSTO GENCODE Reference Quantification
################################################################################

# Import GENCODE v45 transcript-to-gene mapping
tx2g.gencode <- data.frame(import("/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v45/gencode.v45.annotation.gtf"))
tx2g.gencode <- tx2g.gencode[tx2g.gencode$seqnames %in% 
                               c(paste0("chr", seq.int(1, 22, 1)), "chrX", "chrY", "chrM"), ]
tx2g.gencode <- tx2g.gencode[, c("transcript_id", "gene_id")]
tx2g.gencode <- tx2g.gencode[!is.na(tx2g.gencode$transcript_id), ]
tx2g.gencode <- tx2g.gencode[!duplicated(tx2g.gencode), ]

# Process GENCODE quantification data
dirs.salmon.gencode <- list.dirs("/rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/infrv_GUSTO/COUNTS_gencode", recursive = FALSE)
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
# GUSTO GENCODE+ (combined) Reference Quantification
################################################################################

# Import GENCODE+ (combined) v45 transcript-to-gene mapping
tx2g.combined <- data.frame(import("Assembly/lr_assembly_GENCODE_v45_combined.gtf"))
tx2g.combined <- tx2g.combined[tx2g.combined$seqnames %in% 
                                 c(paste0("chr", seq.int(1, 22, 1)), "chrX", "chrY", "chrM"), ]
tx2g.combined <- tx2g.combined[, c("transcript_id", "gene_id")]
tx2g.combined <- tx2g.combined[!is.na(tx2g.combined$transcript_id), ]
tx2g.combined <- tx2g.combined[!duplicated(tx2g.combined), ]

# Process GENCODE+ (combined) quantification data
dirs.salmon.combined <- list.dirs("/rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/infrv_GUSTO/COUNTS_combined", recursive = FALSE)
s.combined <- catchSalmon(dirs.salmon.combined)

files.list.combined <- list()
for (i in 1:length(dirs.salmon.combined)) {
  files.list.combined[i] <- paste(unlist(dirs.salmon.combined[i]), "quant.sf", sep = "/")
}
dirnames.combined <- sapply(strsplit(unlist(dirs.salmon.combined), "/"), "[",
                            length(unlist(strsplit(unlist(dirs.salmon.combined), "/")[1])))
names(files.list.combined) <- dirnames.combined

# Import and process GENCODE+ (combined) quantification with QU correction
txi.combined <- tximport(unlist(files.list.combined), 
                         type = "salmon", 
                         tx2gene = tx2g.combined,
                         geneIdCol = "gene_id", 
                         txIdCol = "transcript_id", 
                         dropInfReps = TRUE, 
                         txOut = TRUE)

txi.combined$counts <- txi.combined$counts / s.combined$annotation$Overdispersion
txi.combined$abundance <- cpm(txi.combined$counts)
txi.combined.GUSTO <- txi.combined

dds.combined.GUSTO <- DESeqDataSetFromTximport(txi = txi.combined.GUSTO, 
                                               colData = cordat.GUSTO, 
                                               design = ~1)

################################################################################
# Save tx2g tables for downstream analyses
################################################################################
save(list=c("tx2g.assembly","tx2g.gencode","tx2g.combined"),
     file="Supplementary/tx2g_ESPRESSO_assembly_SC_filtered.RData")

################################################################################
# Gen3G Cohort Analysis
################################################################################

# Prepare Gen3G covariate data
cordat <- read.csv("/rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/metadata_Gen3G.csv") # From dbGaP
cordat$GA <- as.numeric(scale(cordat$GA))

names(cordat)[which(names(cordat) == "SEX")] <- "sex"
cordat$sex <- factor(cordat$sex)

names(cordat)[which(names(cordat) == "Visite_Diagnostic_de_DG_V2")] <- "gdm"
cordat <- cordat[!cordat$gdm == 999, ]
cordat$gdm <- factor(cordat$gdm)

names(cordat)[which(names(cordat) == "Accouchement_Poids")] <- "birth_weight"
cordat$birth_weight <- scale(cordat$birth_weight)

cordat.Gen3G <- cordat
save(list=c("cordat.Gen3G"),file="Supplementary/metadata_Gen3G.RData")

################################################################################
# Gen3G Assembly Analysis
################################################################################

## Prepare assembly gene expression data
dirs.salmon.assembly <- list.dirs("/rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/infrv_Gen3G/COUNTS_assembly_SC", 
                                  recursive = FALSE)
files.list.assembly <- list()
for (i in 1:length(dirs.salmon.assembly)) {
  files.list.assembly[i] <- paste(unlist(dirs.salmon.assembly[i]), "quant.sf", sep = "/")
}
dirnames.assembly <- sapply(strsplit(unlist(dirs.salmon.assembly), "/"), "[",
                            length(unlist(strsplit(unlist(dirs.salmon.assembly), "/")[1])))
names(files.list.assembly) <- dirnames.assembly

# Filter to samples present in metadata
keep <- which(dirnames.assembly %in% cordat.Gen3G$Run)
dirs.salmon.assembly <- dirs.salmon.assembly[keep]
dirnames.assembly <- dirnames.assembly[keep]
files.list.assembly <- unlist(files.list.assembly)[keep]

# Import gene-level quantification
txi.assembly.g <- tximport(unlist(files.list.assembly), type = "salmon", tx2gene = tx2g.assembly,
                           geneIdCol = "gene_id", txIdCol = "transcript_id", dropInfReps = TRUE, txOut = FALSE)

### Quality control and batch filtering
dds.assembly <- DESeqDataSetFromTximport(txi.assembly.g, cordat.Gen3G, ~1)
vst.assembly <- vst(dds.assembly, blind = TRUE)
pca.assembly.check <- plotPCA(vst.assembly, intgroup = "sex")
pca.assembly.check

# Check batch effects
boxplot(gdm ~ sequencing_batch, cordat.Gen3G)
data.frame(table(cordat.Gen3G[, c("gdm", "sequencing_batch")]))

# Filter out problematic sequencing batches
filter.assembly <- cordat.Gen3G[!cordat.Gen3G$sequencing_batch %in% c("SK-3KKV", "SK-3KKW"), "Run"]
txi.assembly.filter <- filter_tximport(txi.assembly.g, filter.assembly)
cordat.Gen3G <- cordat.Gen3G[cordat.Gen3G$sequencing_batch %in% c("SK-3KKV", "SK-3KKW"), ]

# Save metadata for later use
metadata.Gen3G <- metadata
cordat.Gen3G.Gen3G <- cordat.Gen3G
save(list=c("metadata.Gen3G","cordat.Gen3G.Gen3G"),file="Supplementary/metadata_Gen3G.RData")

################################################################################
# Gen3G Quantification Processing (Assembly, GENCODE, and GENCODE+)
################################################################################

## Gen3G Assembly quantification
dirs.salmon.assembly <- list.dirs("/rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/infrv_Gen3G/COUNTS_assembly_SC", 
                                  recursive = FALSE)
names(dirs.salmon.assembly) <- basename(dirs.salmon.assembly)
dirs.salmon.assembly <- dirs.salmon.assembly[names(dirs.salmon.assembly) %in% cordat.Gen3G$Run]

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
                                               colData = cordat.Gen3G, 
                                               design = ~1)

################################################################################
# Figure S6A - Short-read sample PCAs
################################################################################

# GUSTO
vst.GUSTO <- assay(varianceStabilizingTransformation(dds.assembly.GUSTO))
pca.GUSTO <- prcomp(t(vst.GUSTO))
pca.df.GUSTO <- data.frame(pca.GUSTO$x[, 1:2])
pca.df.GUSTO$SampleID <- row.names(pca.df.GUSTO)
pca.df.GUSTO <- left_join(pca.df.GUSTO, cordat.GUSTO[, c("SampleID", "gdm")])
pca.df.GUSTO[is.na(pca.df.GUSTO)] <- 0
var_explained <- pca.GUSTO$sdev^2
pve <- var_explained / sum(var_explained)
pve[1:2] * 100  # PC1: 25.72%, PC2: 18.26%

p_gusto <- ggplot(
  subset(pca.df, Cohort == "GUSTO (n=200)"),
  aes(x = PC1, y = PC2, color = gdm)
) +
  geom_point() +
  theme_minimal() +
  ggtitle("GUSTO (n=200)") +
  xlab("PC1: 25.72% variance") +
  ylab("PC2: 18.26% variance") +
  theme(legend.position = "none") 

# Gen3G
vst.Gen3G <- assay(varianceStabilizingTransformation(dds.assembly.Gen3G))
pca.Gen3G <- prcomp(t(vst.Gen3G))
pca.df.Gen3G <- data.frame(pca.Gen3G$x[, 1:2])
pca.df.Gen3G$Run <- row.names(pca.df.Gen3G)
pca.df.Gen3G <- left_join(pca.df.Gen3G, cordat.Gen3G[, c("Run", "gdm")])
pca.df.Gen3G[is.na(pca.df.Gen3G)] <- 0
var_explained <- pca.Gen3G$sdev^2
pve <- var_explained / sum(var_explained)
pve[1:2] * 100  # PC1: 19.54%, PC2: 6.36%

p_gen3g <- ggplot(
  subset(pca.df, Cohort == "Gen3G (n=152)"),
  aes(x = PC1, y = PC2, color = gdm)
) +
  geom_point() +
  theme_minimal() +
  ggtitle("Gen3G (n=152)") +
  xlab("PC1: 19.54% variance") +
  ylab("PC2: 6.36% variance") +
  theme(legend.position = "none") 

# Combine
figs6a <- (p_gusto | p_gen3g) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
figs6a <- figs6a + plot_annotation(
    title = "Short-read samples"
  )
figs6a

ggsave("Figures/figs6a.pdf", figs6a, width = 4.15, height = 3.5, units = "in", dpi = 300)

## Gen3G GENCODE quantification
dirs.salmon.gencode <- list.dirs("/rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/infrv_Gen3G/COUNTS_gencode", 
                                 recursive = FALSE)
names(dirs.salmon.gencode) <- basename(dirs.salmon.gencode)
dirs.salmon.gencode <- dirs.salmon.gencode[names(dirs.salmon.gencode) %in% cordat.Gen3G$Run]

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
                                              colData = cordat.Gen3G, 
                                              design = ~1)

## Gen3G GENCODE+ (combined) quantification
dirs.salmon.combined <- list.dirs("/rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/infrv_Gen3G/COUNTS_combined", 
                                  recursive = FALSE)
names(dirs.salmon.combined) <- basename(dirs.salmon.combined)
dirs.salmon.combined <- dirs.salmon.combined[names(dirs.salmon.combined) %in% cordat.Gen3G$Run]

s.combined <- catchSalmon(dirs.salmon.combined)

files.list.combined <- list()
for (i in 1:length(dirs.salmon.combined)) {
  files.list.combined[i] <- paste(unlist(dirs.salmon.combined[i]), "quant.sf", sep = "/")
}
dirnames.combined <- sapply(strsplit(unlist(dirs.salmon.combined), "/"), "[",
                            length(unlist(strsplit(unlist(dirs.salmon.combined), "/")[1])))
names(files.list.combined) <- dirnames.combined

txi.combined <- tximport(unlist(files.list.combined), 
                         type = "salmon", 
                         tx2gene = tx2g.combined,
                         geneIdCol = "gene_id", 
                         txIdCol = "transcript_id", 
                         dropInfReps = TRUE, 
                         txOut = TRUE)

txi.combined$counts <- txi.combined$counts / s.combined$annotation$Overdispersion
txi.combined$abundance <- cpm(txi.combined$counts)
txi.combined.Gen3G <- txi.combined
rm(txi.combined)

dds.combined.Gen3G <- DESeqDataSetFromTximport(txi = txi.combined.Gen3G, 
                                               colData = cordat.Gen3G, 
                                               design = ~1)

################################################################################
# Figure 3B: Total Expression Comparison Across Cohorts and Annotations
################################################################################

## Calculate normalized TPM-like expression values for comparison
# GUSTO Assembly
fig2b.dat.1 <- counts(dds.assembly.GUSTO, normalized = FALSE)
fig2b.dat.1 <- fig2b.dat.1 / txi.assembly.GUSTO$length
fig2b.dat.1 <- log1p(t(t(fig2b.dat.1) * 1e6 / colSums(fig2b.dat.1)))
fig2b.dat.1 <- data.frame(colSums(fig2b.dat.1))
fig2b.dat.1$sample_id <- row.names(fig2b.dat.1)
fig2b.dat.1 <- fig2b.dat.1[, c(2, 1)]
names(fig2b.dat.1)[2] <- "assembly"
# GUSTO GENCODE
fig2b.dat.2 <- counts(dds.gencode.GUSTO, normalized = FALSE)
fig2b.dat.2 <- fig2b.dat.2 / txi.gencode.GUSTO$length
fig2b.dat.2 <- log1p(t(t(fig2b.dat.2) * 1e6 / colSums(fig2b.dat.2)))
fig2b.dat.2 <- data.frame(colSums(fig2b.dat.2))
fig2b.dat.2$sample_id <- row.names(fig2b.dat.2)
fig2b.dat.2 <- fig2b.dat.2[, c(2, 1)]
names(fig2b.dat.2)[2] <- "gencode"
# GUSTO Combined (GENCODE+)
fig2b.dat.2b <- counts(dds.combined.GUSTO, normalized = FALSE)
fig2b.dat.2b <- fig2b.dat.2b / txi.combined.GUSTO$length
fig2b.dat.2b <- log1p(t(t(fig2b.dat.2b) * 1e6 / colSums(fig2b.dat.2b)))

### GENCODE_plus_txs <- row.names(fig2b.dat.2b[rowSums(fig2b.dat.2b)>10,])

fig2b.dat.2b <- data.frame(colSums(fig2b.dat.2b))
fig2b.dat.2b$sample_id <- row.names(fig2b.dat.2b)
fig2b.dat.2b <- fig2b.dat.2b[, c(2, 1)]
names(fig2b.dat.2b)[2] <- "combined"
# Combine GUSTO data
fig2b.dat.a <- data.frame(sample_id = fig2b.dat.1$sample_id)
fig2b.dat.a <- left_join(fig2b.dat.a, fig2b.dat.1)
fig2b.dat.a <- left_join(fig2b.dat.a, fig2b.dat.2)
fig2b.dat.a <- left_join(fig2b.dat.a, fig2b.dat.2b)
# Gen3G Assembly
fig2b.dat.3 <- counts(dds.assembly.Gen3G, normalized = FALSE)
fig2b.dat.3 <- fig2b.dat.3 / txi.assembly.Gen3G$length
fig2b.dat.3 <- log1p(t(t(fig2b.dat.3) * 1e6 / colSums(fig2b.dat.3)))
fig2b.dat.3 <- data.frame(colSums(fig2b.dat.3))
fig2b.dat.3$sample_id <- row.names(fig2b.dat.3)
fig2b.dat.3 <- fig2b.dat.3[, c(2, 1)]
names(fig2b.dat.3)[2] <- "assembly"
# Gen3G GENCODE
fig2b.dat.4 <- counts(dds.gencode.Gen3G, normalized = FALSE)
fig2b.dat.4 <- fig2b.dat.4 / txi.gencode.Gen3G$length
fig2b.dat.4 <- log1p(t(t(fig2b.dat.4) * 1e6 / colSums(fig2b.dat.4)))
fig2b.dat.4 <- data.frame(colSums(fig2b.dat.4))
fig2b.dat.4$sample_id <- row.names(fig2b.dat.4)
fig2b.dat.4 <- fig2b.dat.4[, c(2, 1)]
names(fig2b.dat.4)[2] <- "gencode"
# Gen3G Combined (GENCODE+)
fig2b.dat.4b <- counts(dds.combined.Gen3G, normalized = FALSE)
fig2b.dat.4b <- fig2b.dat.4b / txi.combined.Gen3G$length
fig2b.dat.4b <- log1p(t(t(fig2b.dat.4b) * 1e6 / colSums(fig2b.dat.4b)))

### which(row.names(fig2b.dat.4b)%in%GENCODE_plus_txs)
### fig2b.dat.4b.sub <- fig2b.dat.4b[row.names(fig2b.dat.4b)%in%GENCODE_plus_txs,]
# fig2b.dat.4b <- data.frame(colSums(fig2b.dat.4b.sub))

fig2b.dat.4b <- data.frame(colSums(fig2b.dat.4b))
fig2b.dat.4b$sample_id <- row.names(fig2b.dat.4b)
fig2b.dat.4b <- fig2b.dat.4b[, c(2, 1)]
names(fig2b.dat.4b)[2] <- "combined"
# Combine Gen3G data
fig2b.dat.b <- data.frame(sample_id = fig2b.dat.3$sample_id)
fig2b.dat.b <- left_join(fig2b.dat.b, fig2b.dat.3)
fig2b.dat.b <- left_join(fig2b.dat.b, fig2b.dat.4)
fig2b.dat.b <- left_join(fig2b.dat.b, fig2b.dat.4b)

# Create comprehensive comparison dataset
fig2b.dat <- data.frame(
  Cohort = c(rep.int("GUSTO n=200", length(fig2b.dat.a$assembly) * 3),
             rep.int("Gen3G n=152", length(fig2b.dat.b$assembly) * 3)),
  annotation = c(rep.int("lr-assembly", length(fig2b.dat.a$assembly)),
                 rep.int("GENCODEv45", length(fig2b.dat.a$gencode)),
                 rep.int("GENCODE+", length(fig2b.dat.a$combined)),
                 rep.int("lr-assembly", length(fig2b.dat.b$assembly)),
                 rep.int("GENCODEv45", length(fig2b.dat.b$gencode)),
                 rep.int("GENCODE+", length(fig2b.dat.b$combined))),
  sum = c(fig2b.dat.a$assembly, fig2b.dat.a$gencode, fig2b.dat.a$combined,
          fig2b.dat.b$assembly, fig2b.dat.b$gencode, fig2b.dat.b$combined)
)

# Set factor order for annotation
fig2b.dat$annotation <- factor(fig2b.dat$annotation, 
                               levels = c("GENCODEv45", "GENCODE+", "lr-assembly"))
fig2b.dat$Cohort <- sub(" n=([0-9]+)", " (n=\\1)", fig2b.dat$Cohort)

## Generate visualization
fig3b <- ggplot(fig2b.dat, aes(x = Cohort, y = sum, fill = annotation)) +
  geom_violin(alpha = 0.5, position = position_dodge(0.9)) + 
  geom_boxplot(alpha = 0.5, width = 0.1, position = position_dodge(0.9)) + 
  theme_minimal() +
  theme(legend.position = "top",
        axis.title.x = element_blank(),
        legend.title = element_blank()) +
  labs(y = expression(Sigma*"  log(TPM)")) +
  ggtitle("Total isoform expression across samples") +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values = c("royalblue", "#de77ae", "seagreen"))
fig3b

ggsave("Figures/fig3b.pdf", fig3b, width = 4.6, height = 2.46, units = "in", dpi = 300)

## Statistical analysis of expression differences
fig2b.dat$sample <- c(paste0("GUSTO_", 1:200), paste0("GUSTO_", 1:200), paste0("GUSTO_", 1:200),
                      paste0("Gen3G_", 1:152), paste0("Gen3G_", 1:152), paste0("Gen3G_", 1:152))
# Mixed-effects model accounting for paired samples
model <- lmerTest::lmer(sum ~ annotation * Cohort + (1 | sample), data = fig2b.dat)
summary(model)

library(emmeans)
emm <- emmeans(model, ~ annotation | Cohort)
emm

################################################################################
# Expression Pattern Analysis: Mean and Variance by Transcript Categories
################################################################################

## Define transcript categories for comparison
txs.all <- union(tx2g.assembly$transcript_id, tx2g.gencode$transcript_id)
txs.both <- intersect(tx2g.assembly$transcript_id, tx2g.gencode$transcript_id)

# Load raw assembly data to identify filtering effects
gtf.assembly_all <- "Assembly/ESPRESSO_corrected_corrected.gtf"
all_assembled_tx <- unique(data.frame(import(gtf.assembly_all))$transcript_id)

# Classify transcripts by assembly/filtering status
not_assembled <- setdiff(txs.all, all_assembled_tx)
filtered <- setdiff(all_assembled_tx, tx2g.assembly$transcript_id)
shared <- txs.both
novel <- setdiff(tx2g.assembly$transcript_id, txs.both)

################################################################################
# Figure 3C: Mean Expression Distribution by Transcript Categories (GUSTO)
################################################################################

## Calculate mean expression for each transcript across GUSTO samples
fig2c1.dat.1 <- counts(estimateSizeFactors(dds.assembly.GUSTO), normalized = TRUE)
fig2c1.dat.1 <- data.frame(rowMeans(fig2c1.dat.1))
fig2c1.dat.1$transcript_id <- row.names(fig2c1.dat.1)
fig2c1.dat.1 <- fig2c1.dat.1[, c(2, 1)]
names(fig2c1.dat.1)[2] <- "lr-assembly"
fig2c1.dat.2 <- counts(estimateSizeFactors(dds.gencode.GUSTO), normalized = TRUE)
fig2c1.dat.2 <- data.frame(rowMeans(fig2c1.dat.2))
fig2c1.dat.2$transcript_id <- row.names(fig2c1.dat.2)
fig2c1.dat.2 <- fig2c1.dat.2[, c(2, 1)]
names(fig2c1.dat.2)[2] <- "GENCODEv45"
fig2c1.dat.3 <- counts(estimateSizeFactors(dds.combined.GUSTO), normalized = TRUE)
fig2c1.dat.3 <- data.frame(rowMeans(fig2c1.dat.3))
fig2c1.dat.3$transcript_id <- row.names(fig2c1.dat.3)
fig2c1.dat.3 <- fig2c1.dat.3[, c(2, 1)]
names(fig2c1.dat.3)[2] <- "GENCODE+"

# Combine data and classify transcripts
fig2c1.dat.a <- data.frame(transcript_id = txs.all)
fig2c1.dat.a <- left_join(fig2c1.dat.a, fig2c1.dat.1)
fig2c1.dat.a <- left_join(fig2c1.dat.a, fig2c1.dat.2)
fig2c1.dat.a <- left_join(fig2c1.dat.a, fig2c1.dat.3)
fig2c1.dat.a[is.na(fig2c1.dat.a)] <- 0
# Assign category labels
fig2c1.dat.a$set <- NA
fig2c1.dat.a[fig2c1.dat.a$transcript_id %in% shared, "set"] <- "shared"
fig2c1.dat.a[fig2c1.dat.a$transcript_id %in% novel, "set"] <- "novel"
fig2c1.dat.a[fig2c1.dat.a$transcript_id %in% filtered, "set"] <- "filtered"
fig2c1.dat.a[fig2c1.dat.a$transcript_id %in% not_assembled, "set"] <- "absent"
# Convert to long format for visualization
fig2c1.dat.a <- fig2c1.dat.a[, -1] %>%
  pivot_longer(
    cols = c(`lr-assembly`, GENCODEv45, `GENCODE+`),
    names_to = "annotation",
    values_to = "mean"
  )
# Remove impossible combinations (novel transcripts can't be in GENCODE, etc.)
fig2c1.dat.a$key <- paste0(fig2c1.dat.a$set, fig2c1.dat.a$annotation)
fig2c1.dat.a <- fig2c1.dat.a[!fig2c1.dat.a$key == "novelGENCODEv45", ]
fig2c1.dat.a <- fig2c1.dat.a[!fig2c1.dat.a$key == "absentlr-assembly", ]
fig2c1.dat.a <- fig2c1.dat.a[!fig2c1.dat.a$key == "filteredlr-assembly", ]
fig2c1.dat.a$set <- factor(fig2c1.dat.a$set, levels = c("shared", "novel", "filtered", "absent"))
fig2c1.dat.a$annotation <- factor(fig2c1.dat.a$annotation, 
                                  levels = c("GENCODEv45", "GENCODE+", "lr-assembly"))
fig2c1.dat.a <- fig2c1.dat.a[order(fig2c1.dat.a$set), ]
fig2c1.dat.a$key <- NULL
fig2c1.dat.a <- fig2c1.dat.a[complete.cases(fig2c1.dat.a), ]

## Create visualization
fig2c1.a <- ggplot(fig2c1.dat.a,
                   aes(x = set, y = mean + 1, fill = annotation)) +
  geom_boxplot(alpha = 0.5, outliers = FALSE) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(),
        panel.grid.minor.y = element_line(),
        legend.title = element_blank(),
        legend.position = "top",
        axis.title.x=element_blank()) +
  ylab("Mean across samples") +
  scale_fill_manual(values = c("royalblue", "#de77ae", "seagreen")) +
  scale_y_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +  
  ggtitle("Isoform expression in lr-assembly vs GENCODE") +
  labs(subtitle="Data from GUSTO (n=200)")
fig2c1.a

################################################################################
# Figure S6B: Mean Expression Distribution (Gen3G)
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
figs4ba.dat.5 <- counts(estimateSizeFactors(dds.combined.Gen3G), normalized = TRUE)
figs4ba.dat.5 <- data.frame(rowMeans(figs4ba.dat.5))
figs4ba.dat.5$transcript_id <- row.names(figs4ba.dat.5)
figs4ba.dat.5 <- figs4ba.dat.5[, c(2, 1)]
names(figs4ba.dat.5)[2] <- "GENCODE+"
# Process Gen3G data with same classification scheme
figs4ba.dat.b <- data.frame(transcript_id = txs.all)
figs4ba.dat.b <- left_join(figs4ba.dat.b, figs4ba.dat.3)
figs4ba.dat.b <- left_join(figs4ba.dat.b, figs4ba.dat.4)
figs4ba.dat.b <- left_join(figs4ba.dat.b, figs4ba.dat.5)
figs4ba.dat.b[is.na(figs4ba.dat.b)] <- 0
figs4ba.dat.b$set <- NA
figs4ba.dat.b[figs4ba.dat.b$transcript_id %in% shared, "set"] <- "shared"
figs4ba.dat.b[figs4ba.dat.b$transcript_id %in% novel, "set"] <- "novel"
figs4ba.dat.b[figs4ba.dat.b$transcript_id %in% filtered, "set"] <- "filtered"
figs4ba.dat.b[figs4ba.dat.b$transcript_id %in% not_assembled, "set"] <- "absent"
figs4ba.dat.b <- figs4ba.dat.b[, -1] %>%
  pivot_longer(
    cols = c(`lr-assembly`, GENCODEv45, `GENCODE+`),
    names_to = "annotation",
    values_to = "mean"
  )
figs4ba.dat.b$key <- paste0(figs4ba.dat.b$set, figs4ba.dat.b$annotation)
figs4ba.dat.b <- figs4ba.dat.b[!figs4ba.dat.b$key == "novelGENCODEv45", ]
figs4ba.dat.b <- figs4ba.dat.b[!figs4ba.dat.b$key == "absentlr-assembly", ]
figs4ba.dat.b <- figs4ba.dat.b[!figs4ba.dat.b$key == "filteredlr-assembly", ]
figs4ba.dat.b$set <- factor(figs4ba.dat.b$set, levels = c("shared", "novel", "filtered", "absent"))
figs4ba.dat.b$annotation <- factor(figs4ba.dat.b$annotation, 
                                   levels = c("GENCODEv45", "GENCODE+", "lr-assembly"))
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
        legend.position = "top",
        axis.title.x=element_blank()) +
  ylab("Mean across samples") +
  scale_fill_manual(values = c("royalblue", "#de77ae", "seagreen")) +
  scale_y_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +  
  ggtitle("Isoform expression in lr-assembly vs GENCODE") +
  labs(subtitle="Data from Gen3G (n=152)")
figs4ba.b

################################################################################
# Statistical Analysis of Mean Expression Differences
################################################################################

## Combine both cohorts for statistical testing
fig2c1.dat.a <- data.frame(fig2c1.dat.a)
fig2c1.dat.a$cohort <- "GUSTO"
figs4ba.dat.b <- data.frame(figs4ba.dat.b)
figs4ba.dat.b$cohort <- "Gen3G"
fig2c1.dat.all <- rbind(fig2c1.dat.a, figs4ba.dat.b)
# Prepare factors for modeling
fig2c1.dat.all$set <- factor(fig2c1.dat.all$set, levels = c("shared", "novel", "filtered", "absent"))
fig2c1.dat.all$annotation <- factor(fig2c1.dat.all$annotation, 
                                    levels = c("GENCODEv45", "GENCODE+", "lr-assembly"))
fig2c1.dat.all$cohort <- factor(fig2c1.dat.all$cohort)
# Nested linear model to test differences between categories
fit_nested <- lm(log1p(mean) ~ cohort:set + cohort:annotation, data = fig2c1.dat.all)
anova(fit_nested)
# Post-hoc pairwise comparisons
emm_nested <- emmeans(fit_nested, ~ set | cohort)
pairwise_contrasts <- contrast(emm_nested, method = "pairwise", adjust = "tukey")
summary(pairwise_contrasts, infer = c(TRUE, TRUE))

################################################################################
# Expression Variance Analysis
################################################################################

## Figure 1H: Variance Distribution by Transcript Categories (GUSTO)
fig2c2.dat.1 <- counts(estimateSizeFactors(dds.assembly.GUSTO), normalized = TRUE)
fig2c2.dat.1 <- data.frame(rowVars(fig2c2.dat.1))
fig2c2.dat.1$transcript_id <- row.names(fig2c2.dat.1)
fig2c2.dat.1 <- fig2c2.dat.1[, c(2, 1)]
names(fig2c2.dat.1)[2] <- "lr-assembly"

fig2c2.dat.2 <- counts(estimateSizeFactors(dds.gencode.GUSTO), normalized = TRUE)
fig2c2.dat.2 <- data.frame(rowVars(fig2c2.dat.2))
fig2c2.dat.2$transcript_id <- row.names(fig2c2.dat.2)
fig2c2.dat.2 <- fig2c2.dat.2[, c(2, 1)]
names(fig2c2.dat.2)[2] <- "GENCODEv45"

fig2c2.dat.3 <- counts(estimateSizeFactors(dds.combined.GUSTO), normalized = TRUE)
fig2c2.dat.3 <- data.frame(rowVars(fig2c2.dat.3))
fig2c2.dat.3$transcript_id <- row.names(fig2c2.dat.3)
fig2c2.dat.3 <- fig2c2.dat.3[, c(2, 1)]
names(fig2c2.dat.3)[2] <- "GENCODE+"

# Process variance data with same classification approach
fig2c2.dat.a <- data.frame(transcript_id = txs.all)
fig2c2.dat.a <- left_join(fig2c2.dat.a, fig2c2.dat.1)
fig2c2.dat.a <- left_join(fig2c2.dat.a, fig2c2.dat.2)
fig2c2.dat.a <- left_join(fig2c2.dat.a, fig2c2.dat.3)
fig2c2.dat.a[is.na(fig2c2.dat.a)] <- 0

fig2c2.dat.a$set <- NA
fig2c2.dat.a[fig2c2.dat.a$transcript_id %in% shared, "set"] <- "shared"
fig2c2.dat.a[fig2c2.dat.a$transcript_id %in% novel, "set"] <- "novel"
fig2c2.dat.a[fig2c2.dat.a$transcript_id %in% filtered, "set"] <- "filtered"
fig2c2.dat.a[fig2c2.dat.a$transcript_id %in% not_assembled, "set"] <- "absent"

fig2c2.dat.a <- fig2c2.dat.a[, -1] %>%
  pivot_longer(
    cols = c(`lr-assembly`, GENCODEv45, `GENCODE+`),
    names_to = "annotation",
    values_to = "mean"
  )

fig2c2.dat.a$key <- paste0(fig2c2.dat.a$set, fig2c2.dat.a$annotation)
fig2c2.dat.a <- fig2c2.dat.a[!fig2c2.dat.a$key == "novelGENCODEv45", ]
fig2c2.dat.a <- fig2c2.dat.a[!fig2c2.dat.a$key == "absentlr-assembly", ]
fig2c2.dat.a <- fig2c2.dat.a[!fig2c2.dat.a$key == "filteredlr-assembly", ]

fig2c2.dat.a$set <- factor(fig2c2.dat.a$set, levels = c("shared", "novel", "filtered", "absent"))
fig2c2.dat.a$annotation <- factor(fig2c2.dat.a$annotation, 
                                  levels = c("GENCODEv45", "GENCODE+", "lr-assembly"))
fig2c2.dat.a <- fig2c2.dat.a[order(fig2c2.dat.a$set), ]
fig2c2.dat.a$key <- NULL
fig2c2.dat.a <- fig2c2.dat.a[complete.cases(fig2c2.dat.a), ]

fig2c2.a <- ggplot(fig2c2.dat.a,
                   aes(x = set, y = mean + 1, fill = annotation)) +
  geom_boxplot(alpha = 0.5, outliers = FALSE) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(),
        panel.grid.minor.y = element_line(),
        legend.title = element_blank(),
        legend.position = "top",
        axis.title.x=element_blank()) +
  ylab("Variance across samples") +
  scale_fill_manual(values = c("royalblue", "#de77ae", "seagreen")) +
  scale_y_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +  
  ggtitle("Isoform expression in lr-assembly vs GENCODE") +
  labs(subtitle="Data from GUSTO (n=200)")

fig2c2.a

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

figs4bb.dat.5 <- counts(estimateSizeFactors(dds.combined.Gen3G), normalized = TRUE)
figs4bb.dat.5 <- data.frame(rowVars(figs4bb.dat.5))
figs4bb.dat.5$transcript_id <- row.names(figs4bb.dat.5)
figs4bb.dat.5 <- figs4bb.dat.5[, c(2, 1)]
names(figs4bb.dat.5)[2] <- "GENCODE+"

figs4bb.dat.b <- data.frame(transcript_id = txs.all)
figs4bb.dat.b <- left_join(figs4bb.dat.b, figs4bb.dat.3)
figs4bb.dat.b <- left_join(figs4bb.dat.b, figs4bb.dat.4)
figs4bb.dat.b <- left_join(figs4bb.dat.b, figs4bb.dat.5)
figs4bb.dat.b[is.na(figs4bb.dat.b)] <- 0

figs4bb.dat.b$set <- NA
figs4bb.dat.b[figs4bb.dat.b$transcript_id %in% shared, "set"] <- "shared"
figs4bb.dat.b[figs4bb.dat.b$transcript_id %in% novel, "set"] <- "novel"
figs4bb.dat.b[figs4bb.dat.b$transcript_id %in% filtered, "set"] <- "filtered"
figs4bb.dat.b[figs4bb.dat.b$transcript_id %in% not_assembled, "set"] <- "absent"

figs4bb.dat.b <- figs4bb.dat.b[, -1] %>%
  pivot_longer(
    cols = c(`lr-assembly`, GENCODEv45, `GENCODE+`),
    names_to = "annotation",
    values_to = "mean"
  )

figs4bb.dat.b$key <- paste0(figs4bb.dat.b$set, figs4bb.dat.b$annotation)
figs4bb.dat.b <- figs4bb.dat.b[!figs4bb.dat.b$key == "novelGENCODEv45", ]
figs4bb.dat.b <- figs4bb.dat.b[!figs4bb.dat.b$key == "absentlr-assembly", ]
figs4bb.dat.b <- figs4bb.dat.b[!figs4bb.dat.b$key == "filteredlr-assembly", ]

figs4bb.dat.b$set <- factor(figs4bb.dat.b$set, levels = c("shared", "novel", "filtered", "absent"))
figs4bb.dat.b$annotation <- factor(figs4bb.dat.b$annotation, 
                                   levels = c("GENCODEv45", "GENCODE+", "lr-assembly"))
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
        legend.position = "top",
        axis.title.x=element_blank()) +
  ylab("Variance across samples") +
  scale_fill_manual(values = c("royalblue", "#de77ae", "seagreen")) +
  scale_y_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +  
  ggtitle("Isoform expression in lr-assembly vs GENCODE") +
  labs(subtitle="Data from Gen3G (n=152)")

figs4bb.b

## Statistical analysis of variance patterns
fig2c2 <- data.frame(fig2c2.dat.a)
fig2c2$cohort <- "GUSTO"
figs4bb.dat.b <- data.frame(figs4bb.dat.b)
figs4bb.dat.b$cohort <- "Gen3G"
fig2c2all <- rbind(fig2c2, figs4bb.dat.b)

fig2c2all$set <- factor(fig2c2all$set, levels = c("shared", "novel", "filtered", "absent"))
fig2c2all$annotation <- factor(fig2c2all$annotation, 
                               levels = c("GENCODEv45", "GENCODE+", "lr-assembly"))
fig2c2all$cohort <- factor(fig2c2all$cohort)
names(fig2c2all)[3] <- "variance"

fit_nested <- lm(log1p(variance) ~ cohort:set + cohort:annotation, data = fig2c2all)
anova(fit_nested)

emm_nested <- emmeans(fit_nested, ~ set | cohort)
pairwise_contrasts <- contrast(emm_nested, method = "pairwise", adjust = "tukey")
summary(pairwise_contrasts, infer = c(TRUE, TRUE))

################################################################################
# Combine panels for fig3c and figs6b
################################################################################

fig2c1.dat.a$metric <- "Mean"
fig2c2.dat.a$metric <- "Variance"
fig2c.dat <- dplyr::bind_rows(fig2c1.dat.a, fig2c2.dat.a)
fig2c.dat$group <- ifelse(fig2c.dat$set %in% c("shared", "novel"),
                          "Discovery",
                          "Filtration")
fig2c.dat$metric <- factor(fig2c.dat$metric, levels = c("Mean", "Variance"))
fig2c.dat$group  <- factor(fig2c.dat$group, levels = c("Discovery", "Filtration"))

fig3c <- ggplot(fig2c.dat,
                aes(x = set, y = mean + 1, fill = annotation)) +
  geom_boxplot(alpha = 0.5, outliers = FALSE) +
  facet_grid(metric ~ group, scales = "free", space = "free_x") +
  theme_classic() +
  theme(
    panel.grid.major.y = element_line(),
    panel.grid.minor.y = element_line(),
    legend.title = element_blank(),
    legend.position = "top",
    axis.title.x = element_blank(),
    strip.background = element_rect(color = "black", fill = NA),
    strip.text = element_text(color = "black")
  ) +
  scale_fill_manual(values = c("royalblue", "#de77ae", "seagreen")) +
  scale_y_log10(labels = scales::label_number(
    scale_cut = scales::cut_short_scale()
  )) +
  ggtitle("Isoform expression changes with transcript discovery and filtration") +
  labs(subtitle = "Data from GUSTO (n=200)", y = NULL)
fig3c

ggsave("Figures/fig3c.pdf", fig3c, width = 5.76, height = 3.68, units = "in", dpi = 300)


figs4ba.dat.b$metric <- "Mean"
figs4bb.dat.b$metric <- "Variance"
figs4.dat <- dplyr::bind_rows(figs4ba.dat.b, figs4bb.dat.b)
figs4.dat$group <- ifelse(figs4.dat$set %in% c("shared", "novel"),
                          "Discovery",
                          "Filtration")
figs4.dat$metric <- factor(figs4.dat$metric, levels = c("Mean", "Variance"))
figs4.dat$group  <- factor(figs4.dat$group, levels = c("Discovery", "Filtration"))

figs6b <- ggplot(figs4.dat,
                aes(x = set, y = mean + 1, fill = annotation)) +
  geom_boxplot(alpha = 0.5, outliers = FALSE) +
  facet_grid(metric ~ group, scales = "free", space = "free_x") +
  theme_classic() +
  theme(
    panel.grid.major.y = element_line(),
    panel.grid.minor.y = element_line(),
    legend.title = element_blank(),
    legend.position = "top",
    axis.title.x = element_blank(),
    strip.background = element_rect(color = "black", fill = NA),
    strip.text = element_text(color = "black")
  ) +
  scale_fill_manual(values = c("royalblue", "#de77ae", "seagreen")) +
  scale_y_log10(labels = scales::label_number(
    scale_cut = scales::cut_short_scale()
  )) +
  ggtitle("Isoform expression changes with transcript discovery and filtration") +
  labs(subtitle = "Data from Gen3G (n=152)", y = NULL)
figs6b

ggsave("Figures/figs6b.pdf", figs6b, width = 5.76, height = 3.68, units = "in", dpi = 300)

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

DIR_COUNTS_assembly="/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/infrv_GUSTO/COUNTS_assembly_SC"
DIR_COUNTS_gencode="/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/infrv_GUSTO/COUNTS_gencode"
DIR_COUNTS_combined="/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/infrv_GUSTO/COUNTS_combined"
DIR_COUNTS_adipose="/rsrch5/home/epi/stbresnahan/scratch/natcomms_revision/COUNTS_GUSTO_Adipose_Subcutaneous"
DIR_COUNTS_fibroblasts="/rsrch5/home/epi/stbresnahan/scratch/natcomms_revision/COUNTS_GUSTO_Cells_Cultured_fibroblasts"

# Prepare file paths for tximeta
coldata.x <- cordat.GUSTO
coldata.x$files <- paste0(DIR_COUNTS_assembly, "/", coldata.x$SampleID, "/quant.sf")
names(coldata.x)[1] <- "names"

coldata.y <- cordat.GUSTO
coldata.y$files <- paste0(DIR_COUNTS_gencode, "/", coldata.y$SampleID, "/quant.sf")
names(coldata.y)[1] <- "names"

coldata.z <- cordat.GUSTO
coldata.z$files <- paste0(DIR_COUNTS_combined, "/", coldata.z$SampleID, "/quant.sf")
names(coldata.z)[1] <- "names"

coldata.a <- cordat.GUSTO
coldata.a$files <- paste0(DIR_COUNTS_adipose, "/", coldata.a$SampleID, "/quant.sf")
names(coldata.a)[1] <- "names"

coldata.f <- cordat.GUSTO
coldata.f$files <- paste0(DIR_COUNTS_fibroblasts, "/", coldata.f$SampleID, "/quant.sf")
names(coldata.f)[1] <- "names"

## Assembly inferential variance calculation
setTximetaBFC("/rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/SQANTI3/ESPRESSO_corrected_SC_filtered")
makeLinkedTxome(indexDir = "/rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/SQANTI3/ESPRESSO_corrected_SC_filtered/salmon_index/ESPRESSO_corrected_SC_filtered",
                source = "ESPRESSO assembly of placenta lr-rnaseq (filtered)",
                organism = "Homo sapiens",
                genome = "GRCh38",
                release = "p14",
                fasta = "/rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/SQANTI3/ESPRESSO_corrected_SC_filtered/salmon_index/gentrome.fa",
                gtf = "/rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/SQANTI3/ESPRESSO_corrected_SC_filtered/ESPRESSO_corrected_SC_filtered.gtf",
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
setTximetaBFC("/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode.v45.salmon_index")
makeLinkedTxome(indexDir = "/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode.v45.salmon_index/gencode_v45",
                source = "GENCODEv45",
                organism = "Homo sapiens",
                genome = "GRCh38",
                release = "p14",
                fasta = "/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode.v45.salmon_index/gencode.v45.salmon_index/gentrome.fa",
                gtf = "/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v45/gencode.v45.annotation.gtf",
                write = FALSE)

se.y <- tximeta(coldata.y[, c(53, 1:52)], type = "salmon")
se.y <- labelKeep(se.y)
se.y <- se.y[mcols(se.y)$keep, ]

infrv.gencode.GUSTO <- computeInfRV(se.y)
infrv.gencode.GUSTO <- mcols(infrv.gencode.GUSTO)$meanInfRV
rm(se.y)
gc()

## GENCODE+ inferential variance calculation
setTximetaBFC("/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/combined_GTF/salmon")
makeLinkedTxome(indexDir = "/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/combined_GTF/salmon/combined",
                source = "GENCODE+",
                organism = "Homo sapiens",
                genome = "GRCh38",
                release = "p14",
                fasta = "/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/combined_GTF/GCA_000001405.15_GRCh38_no_alt_analysis_set_cleaned_ready_for_salmon_combined.fasta",
                gtf = "/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/combined_GTF/lr_assembly_GENCODE_v45_combined.gtf",
                write = FALSE)

se.z <- tximeta(coldata.z[, c(53, 1:52)], type = "salmon")
se.z <- labelKeep(se.z)
se.z <- se.z[mcols(se.z)$keep, ]

infrv.combined.GUSTO <- computeInfRV(se.z)
infrv.combined.GUSTO <- mcols(infrv.combined.GUSTO)$meanInfRV
rm(se.z)
gc()

## Adipose inferential variance calculation
setTximetaBFC("/rsrch5/home/epi/stbresnahan/scratch/natcomms_revision/tximeta")
makeLinkedTxome(indexDir = "/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/GTEx_v9/SQANTI3/Adipose_Subcutaneous_filtered/salmon/Adipose_Subcutaneous",
                source = "Adipose",
                organism = "Homo sapiens",
                genome = "GRCh38",
                release = "p14",
                fasta = "/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/GTEx_v9/SQANTI3/Adipose_Subcutaneous_filtered/Adipose_Subcutaneous_filtered.fasta",
                gtf = "/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/GTEx_v9/SQANTI3/Adipose_Subcutaneous_filtered/Adipose_Subcutaneous_filtered.fixed.gtf",
                write = FALSE)

se.a <- tximeta(coldata.a[, c(53, 1:52)], type = "salmon")
se.a <- labelKeep(se.a)
se.a <- se.a[mcols(se.a)$keep, ]

infrv.adipose.GUSTO <- computeInfRV(se.a)
infrv.adipose.GUSTO <- mcols(infrv.adipose.GUSTO)$meanInfRV
rm(se.a)
gc()

## Fibroblasts inferential variance calculation
setTximetaBFC("/rsrch5/home/epi/stbresnahan/scratch/natcomms_revision/tximeta")
makeLinkedTxome(indexDir = "/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/GTEx_v9/SQANTI3/Cells_Cultured_fibroblasts_filtered/salmon/Cells_Cultured_fibroblasts",
                source = "Adipose",
                organism = "Homo sapiens",
                genome = "GRCh38",
                release = "p14",
                fasta = "/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/GTEx_v9/SQANTI3/Cells_Cultured_fibroblasts_filtered/Cells_Cultured_fibroblasts_filtered.fasta",
                gtf = "/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/GTEx_v9/SQANTI3/Cells_Cultured_fibroblasts_filtered/Cells_Cultured_fibroblasts_filtered.fixed.gtf",
                write = FALSE)

se.f <- tximeta(coldata.f[, c(53, 1:52)], type = "salmon")
se.f <- labelKeep(se.f)
se.f <- se.f[mcols(se.f)$keep, ]

infrv.fibroblasts.GUSTO <- computeInfRV(se.f)
infrv.fibroblasts.GUSTO <- mcols(infrv.fibroblasts.GUSTO)$meanInfRV
rm(se.f)
gc()

## Process and visualize GUSTO inferential uncertainty
### Union
infrv.GUSTO <- data.frame(transcript_id = union(union(union(union(names(infrv.assembly.GUSTO),
                                                                  names(infrv.gencode.GUSTO)),
                                                            names(infrv.combined.GUSTO)),
                                                      names(infrv.fibroblasts.GUSTO)),
                                                names(infrv.adipose.GUSTO)))

saveRDS(infrv.GUSTO,"infrv_GUSTO.RDS")
# infrv.GUSTO <- readRDS("infrv_GUSTO.RDS")

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

infrv.combined.GUSTO <- data.frame(infrv.combined.GUSTO)
infrv.combined.GUSTO$transcript_id <- row.names(infrv.combined.GUSTO)
infrv.combined.GUSTO <- infrv.combined.GUSTO[, c(2, 1)]
names(infrv.combined.GUSTO)[2] <- "GENCODE+"
infrv.GUSTO <- left_join(infrv.GUSTO, infrv.combined.GUSTO)

infrv.fibroblasts.GUSTO <- data.frame(infrv.fibroblasts.GUSTO)
infrv.fibroblasts.GUSTO$transcript_id <- row.names(infrv.fibroblasts.GUSTO)
infrv.fibroblasts.GUSTO <- infrv.fibroblasts.GUSTO[, c(2, 1)]
names(infrv.fibroblasts.GUSTO)[2] <- "Cells - Cultured Fibroblasts"
infrv.GUSTO <- left_join(infrv.GUSTO, infrv.fibroblasts.GUSTO)

infrv.adipose.GUSTO <- data.frame(infrv.adipose.GUSTO)
infrv.adipose.GUSTO$transcript_id <- row.names(infrv.adipose.GUSTO)
infrv.adipose.GUSTO <- infrv.adipose.GUSTO[, c(2, 1)]
names(infrv.adipose.GUSTO)[2] <- "Adipose - Subcutaneous"
infrv.GUSTO <- left_join(infrv.GUSTO, infrv.adipose.GUSTO)

infrv.GUSTO <- infrv.GUSTO %>%
  pivot_longer(
    cols = c(`lr-assembly`, GENCODEv45, `GENCODE+`, 
             `Cells - Cultured Fibroblasts`, `Adipose - Subcutaneous`),
    names_to = "annotation",
    values_to = "InfRV"
  )

infrv.GUSTO <- infrv.GUSTO[complete.cases(infrv.GUSTO), ]

infrv.GUSTO[infrv.GUSTO$annotation=="lr-assembly","annotation"] <- "Placenta (lr-assembly)"
infrv.GUSTO$annotation <- factor(infrv.GUSTO$annotation, 
                                 levels = c("GENCODEv45", "GENCODE+", "Placenta (lr-assembly)",
                                            "Cells - Cultured Fibroblasts", "Adipose - Subcutaneous"),
                                 labels = c("GENCODEv45", "GENCODE+", "Placenta (lr-assembly)",
                                            "Cells - Cultured Fibroblasts", "Adipose - Subcutaneous"))

################################################################################
# Figure S6I: InfRV in GUSTO + GTEx tissues
################################################################################

figs6i <- ggplot(infrv.GUSTO, 
                aes(x = annotation, y = log10(InfRV), fill = annotation)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.3, width = 0.5, position = position_dodge(0.8)) +
  stat_slab(adjust = .5, width = .6, justification = -.5, 
            position = position_dodge(0.8)) +
  stat_pointinterval(aes(x = annotation, y = log10(InfRV), color = annotation),
                     point_interval = mode_hdi,
                     .width = 0,
                     point_size = 4,
                     shape = 18,
                     linewidth = 0,
                     position = position_dodge(0.8),
                     show.legend = FALSE) +
  scale_fill_manual(values = rev(c("#ff68a1", "#e68613", "seagreen", "#de77ae", "royalblue"))) +
  scale_color_manual(values = rev(c("#ff68a1", "#e68613", "seagreen", "#de77ae", "royalblue"))) +
  scale_x_discrete(limits = rev) +
  theme_classic() +
  labs(fill = "") +
  ggtitle("") +
  theme(panel.grid.major.y = element_line(),
        panel.grid.minor.y = element_line(),
        legend.position = "none") +
  scale_y_continuous(expand = c(0.02, 0)) +
  xlab("") + 
  ylab("log10(mean InfRV)") +
  labs(title="Inferential uncertainty in isoform quantification",
       subtitle="Diamond = mode; data from GUSTO (n=200)") +
  coord_flip()
figs6i

ggsave("Figures/figs6i.pdf", figs6i, width = 5.76, height = 3.68, units = "in", dpi = 300)


################################################################################
# Figure 3D: InfRV in GUSTO
################################################################################

fig2d.dat <- infrv.GUSTO %>% subset(annotation %in% c("GENCODEv45", "GENCODE+", "Placenta (lr-assembly)"))
fig2d.dat$annotation <- factor(fig2d.dat$annotation, 
                                 levels = c("GENCODEv45", "GENCODE+", "Placenta (lr-assembly)"),
                                 labels = c("GENCODEv45", "GENCODE+", "lr-assembly"))

fig3d <- ggplot(fig2d.dat, 
                 aes(x = annotation, y = log10(InfRV), fill = annotation)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.3, width = 0.5, position = position_dodge(0.8)) +
  stat_slab(adjust = .5, width = .6, justification = -.5, 
            position = position_dodge(0.8)) +
  stat_pointinterval(aes(x = annotation, y = log10(InfRV), color = annotation),
                     point_interval = mode_hdi,
                     .width = 0,
                     point_size = 4,
                     shape = 18,
                     linewidth = 0,
                     position = position_dodge(0.8),
                     show.legend = FALSE) +
  scale_fill_manual(values = rev(c("seagreen", "#de77ae", "royalblue"))) +
  scale_color_manual(values = rev(c("seagreen", "#de77ae", "royalblue"))) +
  scale_x_discrete(limits = rev) +
  theme_classic() +
  labs(fill = "") +
  ggtitle("") +
  theme(panel.grid.major.y = element_line(),
        panel.grid.minor.y = element_line(),
        legend.position = "none") +
  scale_y_continuous(expand = c(0.02, 0)) +
  xlab("") + 
  ylab("log10(mean InfRV)") +
  labs(title="Inferential uncertainty in isoform quantification",
       subtitle="Diamond = mode; data from GUSTO (n=200)") +
  coord_flip()
fig3d

ggsave("Figures/fig3d.pdf", fig3d, width = 5.76, height = 3.68, units = "in", dpi = 300)


################################################################################
# Figure S6D: InfRV in GUSTO by structural category
################################################################################


# By structural category in lr-assembly
# Prepare the data - Placenta annotation only
figsf.dat <- infrv.GUSTO %>% 
  subset(annotation == "Placenta (lr-assembly)")

# Merge with structural category
figsf.dat <- merge(figsf.dat, 
                   classif[,c("transcript_id","structural_category")], 
                   by = "transcript_id")

# Create grouped structural category
figsf.dat$structural_category_grouped <- as.character(figsf.dat$structural_category)
figsf.dat$structural_category_grouped[!figsf.dat$structural_category %in% c("FSM", "ISM", "NIC", "NNC")] <- "Other"

# Set factor levels
figsf.dat$structural_category_grouped <- factor(figsf.dat$structural_category_grouped, 
                                                levels = c("FSM", "ISM", "NIC", "NNC", "Other"))

# Color palette
myPalette <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "black")

# Create plot
figs6d <- ggplot(figsf.dat, 
                aes(x = structural_category_grouped, y = log10(InfRV), fill = structural_category_grouped)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.3, width = 0.5, position = position_dodge(0.8)) +
  stat_slab(adjust = .5, width = .6, justification = -.5, 
            position = position_dodge(0.8)) +
  stat_pointinterval(aes(x = structural_category_grouped, y = log10(InfRV), color = structural_category_grouped),
                     point_interval = mode_hdi,
                     .width = 0,
                     point_size = 4,
                     shape = 18,
                     linewidth = 0,
                     position = position_dodge(0.8),
                     show.legend = FALSE) +
  scale_fill_manual(values = myPalette) +
  scale_color_manual(values = myPalette) +
  scale_x_discrete(limits = rev) +
  theme_classic() +
  labs(fill = "") +
  ggtitle("") +
  theme(panel.grid.major.y = element_line(),
        panel.grid.minor.y = element_line(),
        legend.position = "none") +
  scale_y_continuous(expand = c(0.02, 0)) +
  xlab("") + 
  ylab("log10(mean InfRV)") +
  labs(title = "Inferential uncertainty by structural category",
       subtitle = "Diamond = mode; Placenta lr-assembly, GUSTO (n=200)") +
  coord_flip()

figs6d

ggsave("Figures/figs6d.pdf", figs6d, width = 4.82, height = 3.14, units = "in", dpi = 300)


################################################################################
# Gen3G Inferential Uncertainty Analysis
################################################################################

## Setup Gen3G metadata for inferential variance
DIR_COUNTS_assembly <- "/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/infrv_Gen3G/COUNTS_assembly_SC"
DIR_COUNTS_gencode <- "/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/infrv_Gen3G/COUNTS_gencode"
DIR_COUNTS_combined <- "/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/infrv_Gen3G/COUNTS_combined"

# Prepare file paths for Gen3G inferential variance analysis
coldata.x <- cordat.Gen3G
names(coldata.x)[1] <- "names"
coldata.x$files <- paste0(DIR_COUNTS_assembly, "/", coldata.x$names, "/quant.sf")

coldata.y <- cordat.Gen3G
names(coldata.y)[1] <- "names"
coldata.y$files <- paste0(DIR_COUNTS_gencode, "/", coldata.y$names, "/quant.sf")

coldata.z <- cordat.Gen3G
names(coldata.z)[1] <- "names"
coldata.z$files <- paste0(DIR_COUNTS_combined, "/", coldata.z$names, "/quant.sf")

## Gen3G Assembly inferential variance
setTximetaBFC("/rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/SQANTI3/ESPRESSO_corrected_SC_filtered")

se.x <- tximeta(coldata.x[, c(32, 1:31)], type = "salmon")
se.x <- labelKeep(se.x)
se.x <- se.x[mcols(se.x)$keep, ]

infrv.assembly.Gen3G <- computeInfRV(se.x)
infrv.assembly.Gen3G <- mcols(infrv.assembly.Gen3G)$meanInfRV
rm(se.x)
gc()

## Gen3G GENCODE inferential variance
setTximetaBFC("/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode.v45.salmon_index")

se.y <- tximeta(coldata.y[, c(32, 1:31)], type = "salmon")
se.y <- labelKeep(se.y)
se.y <- se.y[mcols(se.y)$keep, ]

infrv.gencode.Gen3G <- computeInfRV(se.y)
infrv.gencode.Gen3G <- mcols(infrv.gencode.Gen3G)$meanInfRV
rm(se.y)
gc()

## Gen3G GENCODE+ inferential variance
setTximetaBFC("/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/combined_GTF/salmon")

se.z <- tximeta(coldata.z[, c(32, 1:31)], type = "salmon")
se.z <- labelKeep(se.z)
se.z <- se.z[mcols(se.z)$keep, ]

infrv.combined.Gen3G <- computeInfRV(se.z)
infrv.combined.Gen3G <- mcols(infrv.combined.Gen3G)$meanInfRV
rm(se.z)
gc()

## Process and visualize Gen3G inferential uncertainty
infrv.Gen3G <- data.frame(transcript_id = union(union(names(infrv.assembly.Gen3G),
                                                      names(infrv.gencode.Gen3G)),
                                                names(infrv.combined.Gen3G)))

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

infrv.combined.Gen3G <- data.frame(infrv.combined.Gen3G)
infrv.combined.Gen3G$transcript_id <- row.names(infrv.combined.Gen3G)
infrv.combined.Gen3G <- infrv.combined.Gen3G[, c(2, 1)]
names(infrv.combined.Gen3G)[2] <- "GENCODE+"
infrv.Gen3G <- left_join(infrv.Gen3G, infrv.combined.Gen3G)

infrv.Gen3G <- infrv.Gen3G %>%
  pivot_longer(
    cols = c(`lr-assembly`, GENCODEv45, `GENCODE+`),
    names_to = "annotation",
    values_to = "InfRV"
  )

infrv.Gen3G <- infrv.Gen3G[complete.cases(infrv.Gen3G), ]
names(infrv.Gen3G)[1] <- "isoform"
infrv.Gen3G <- left_join(infrv.Gen3G, classif.pcorr[, c("isoform", "structural_category")])
infrv.Gen3G$structural_category <- as.character(infrv.Gen3G$structural_category)
infrv.Gen3G[infrv.Gen3G$annotation == "GENCODEv45", "structural_category"] <- "GENCODE"
infrv.Gen3G[infrv.Gen3G$annotation == "GENCODE+" & is.na(infrv.Gen3G$structural_category), "structural_category"] <- "GENCODE"
infrv.Gen3G$structural_category <- factor(infrv.Gen3G$structural_category,
                                          levels = c("GENCODE", "FSM", "ISM", "NIC", "NNC", "Other"))

# Set annotation factor order
infrv.Gen3G$annotation <- factor(infrv.Gen3G$annotation, 
                                 levels = c("GENCODEv45", "GENCODE+", "lr-assembly"))


################################################################################
# Figure S6C: InfRV in Gen3G
################################################################################

## Figure S6c: Inferential uncertainty distribution (Gen3G)
figs6c <- ggplot(infrv.Gen3G, 
                aes(x = annotation, y = log10(InfRV), fill = annotation)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.3, width = 0.5, position = position_dodge(0.8)) +
  stat_slab(adjust = .5, width = .6, justification = -.5, 
            position = position_dodge(0.8)) +
  stat_pointinterval(aes(x = annotation, y = log10(InfRV), color = annotation),
                     point_interval = mode_hdi,
                     .width = 0,
                     point_size = 4,
                     shape = 18,
                     linewidth = 0,
                     position = position_dodge(0.8),
                     show.legend = FALSE) +
  scale_fill_manual(values = rev(c("seagreen", "#de77ae", "royalblue"))) +
  scale_color_manual(values = rev(c("seagreen", "#de77ae", "royalblue"))) +
  scale_x_discrete(limits = rev) +
  theme_classic() +
  labs(fill = "") +
  ggtitle("") +
  theme(panel.grid.major.y = element_line(),
        panel.grid.minor.y = element_line(),
        legend.position = "none") +
  scale_y_continuous(expand = c(0.02, 0)) +
  xlab("") + 
  ylab("log10(mean InfRV)") +
  labs(title="Inferential uncertainty in isoform quantification",
       subtitle="Diamond = mode; data from Gen3G (n=152)") +
  coord_flip()
figs6c

ggsave("Figures/figs6c.pdf", figs6c, width = 5.76, height = 3.68, units = "in", dpi = 300)


################################################################################
# Figure S6E: InfRV in Gen3G by structural category
################################################################################

# By structural category in lr-assembly for Gen3G
# Prepare the data - Placenta annotation only
figs2g.dat <- infrv.Gen3G %>% 
  subset(annotation == "lr-assembly")

# Create grouped structural category
figs2g.dat$structural_category_grouped <- as.character(figs2g.dat$structural_category)
figs2g.dat$structural_category_grouped[!figs2g.dat$structural_category %in% c("FSM", "ISM", "NIC", "NNC")] <- "Other"

# Set factor levels
figs2g.dat$structural_category_grouped <- factor(figs2g.dat$structural_category_grouped, 
                                                 levels = c("FSM", "ISM", "NIC", "NNC", "Other"))

# Color palette
myPalette <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "black")

# Create plot
figs6e <- ggplot(figs2g.dat, 
                 aes(x = structural_category_grouped, y = log10(InfRV), fill = structural_category_grouped)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.3, width = 0.5, position = position_dodge(0.8)) +
  stat_slab(adjust = .5, width = .6, justification = -.5, 
            position = position_dodge(0.8)) +
  stat_pointinterval(aes(x = structural_category_grouped, y = log10(InfRV), color = structural_category_grouped),
                     point_interval = mode_hdi,
                     .width = 0,
                     point_size = 4,
                     shape = 18,
                     linewidth = 0,
                     position = position_dodge(0.8),
                     show.legend = FALSE) +
  scale_fill_manual(values = myPalette) +
  scale_color_manual(values = myPalette) +
  scale_x_discrete(limits = rev) +
  theme_classic() +
  labs(fill = "") +
  ggtitle("") +
  theme(panel.grid.major.y = element_line(),
        panel.grid.minor.y = element_line(),
        legend.position = "none") +
  scale_y_continuous(expand = c(0.02, 0)) +
  xlab("") + 
  ylab("log10(mean InfRV)") +
  labs(title = "Inferential uncertainty by structural category",
       subtitle = "Diamond = mode; Placenta lr-assembly, Gen3G (n=152)") +
  coord_flip()

figs6e
ggsave("Figures/figs6e.pdf", figs6e, width = 4.82, height = 3.14, units = "in", dpi = 300)


################################################################################
# FigS6H: Exon overlaps by InfRV
################################################################################

## Load GTF files
gtf.assembly <- rtracklayer::import("Assembly/ESPRESSO_corrected_SC_filtered.gtf")
gtf.gencode <- rtracklayer::import("/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v45/gencode.v45.annotation.gtf")
gtf.combined <- rtracklayer::import("/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/combined_GTF/lr_assembly_GENCODE_v45_combined.gtf")

## Filter for exons only
exons.assembly <- gtf.assembly[gtf.assembly$type == "exon"]
exons.gencode <- gtf.gencode[gtf.gencode$type == "exon"]
exons.combined <- gtf.combined[gtf.combined$type == "exon"]

## Function to count overlapping exons for a transcript
count_overlapping_exons <- function(transcript_id, exon_gr, min_overlap = 0.8) {
  # Get exons for this transcript
  tx_exons <- exon_gr[exon_gr$transcript_id == transcript_id]
  
  if (length(tx_exons) == 0) return(0)
  
  # Get all other exons from the same gene
  gene_id <- unique(tx_exons$gene_id)
  if (length(gene_id) == 0) return(0)
  
  gene_exons <- exon_gr[exon_gr$gene_id == gene_id & exon_gr$transcript_id != transcript_id]
  
  if (length(gene_exons) == 0) return(0)
  
  # Find all overlaps at once
  overlaps <- findOverlaps(tx_exons, gene_exons)
  
  if (length(overlaps) == 0) return(0)
  
  # Vectorized overlap calculations
  tx_indices <- queryHits(overlaps)
  gene_indices <- subjectHits(overlaps)
  
  tx_widths <- width(tx_exons)[tx_indices]
  gene_widths <- width(gene_exons)[gene_indices]
  
  # Vectorized overlap width calculation
  overlap_starts <- pmax(start(tx_exons)[tx_indices], start(gene_exons)[gene_indices])
  overlap_ends <- pmin(end(tx_exons)[tx_indices], end(gene_exons)[gene_indices])
  overlap_widths <- overlap_ends - overlap_starts + 1
  
  # Vectorized filtering
  valid <- (overlap_widths / tx_widths >= min_overlap) & 
    (overlap_widths / gene_widths >= min_overlap)
  
  return(sum(valid))
}

# Parallel wrapper
count_overlapping_exons_parallel <- function(transcript_ids, exon_gr, min_overlap = 0.8, ncores = parallel::detectCores() - 1) {
  
  # For large exon_gr objects, pre-split by gene to reduce data transfer
  if (object.size(exon_gr) > 50e6) {  # If larger than 50MB
    message("Large exon object detected - using chunked processing")
    
    # Process in chunks to avoid memory issues
    chunk_size <- ceiling(length(transcript_ids) / ncores)
    chunks <- split(transcript_ids, ceiling(seq_along(transcript_ids) / chunk_size))
    
    pb <- txtProgressBar(min = 0, max = length(chunks), style = 3)
    
    results <- list()
    for (i in seq_along(chunks)) {
      results[[i]] <- sapply(chunks[[i]], count_overlapping_exons, 
                             exon_gr = exon_gr, min_overlap = min_overlap)
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    return(unlist(results))
  }
  
  # Original parallel version for smaller objects
  cl <- parallel::makeCluster(ncores)
  on.exit(parallel::stopCluster(cl))
  
  parallel::clusterEvalQ(cl, {
    .libPaths(c("/home/stbresnahan/R/ubuntu/4.3.1", .libPaths()))
    library(GenomicRanges)
  })
  
  parallel::clusterExport(cl, 
                          varlist = c("count_overlapping_exons", "exon_gr", "min_overlap"), 
                          envir = environment())
  
  results <- pbapply::pblapply(transcript_ids, 
                               count_overlapping_exons,
                               exon_gr = exon_gr,
                               min_overlap = min_overlap,
                               cl = cl)
  
  return(unlist(results))
}

## Get full range transcripts (all InfRV values)
all_transcripts.assembly <- infrv.GUSTO.plot %>%
  filter(annotation == "lr-assembly") %>%
  pull(isoform) %>%
  unique() %>% 
  sample(min(1000, length(.)))

all_transcripts.gencode <- infrv.GUSTO.plot %>%
  filter(annotation == "GENCODEv45") %>%
  pull(isoform) %>%
  unique() %>% 
  sample(min(1000, length(.)))

all_transcripts.combined <- infrv.GUSTO.plot %>%
  filter(annotation == "GENCODE+") %>%
  pull(isoform) %>%
  unique() %>% 
  sample(min(1000, length(.)))

# lr-assembly all
overlap_counts.assembly.all <- count_overlapping_exons_parallel(all_transcripts.assembly, exons.assembly)
overlap_data.assembly.all <- data.frame(
  transcript_id = all_transcripts.assembly,
  overlap_count = overlap_counts.assembly.all,
  annotation = "lr-assembly",
  group = "all"
)
closeAllConnections()

# GENCODE all
overlap_counts.gencode.all <- count_overlapping_exons_parallel(all_transcripts.gencode, exons.gencode)
overlap_data.gencode.all <- data.frame(
  transcript_id = all_transcripts.gencode,
  overlap_count = overlap_counts.gencode.all,
  annotation = "GENCODEv45",
  group = "all"
)
closeAllConnections()

# GENCODE+ all
overlap_counts.combined.all <- count_overlapping_exons_parallel(all_transcripts.combined, exons.combined)
overlap_data.combined.all <- data.frame(
  transcript_id = all_transcripts.combined,
  overlap_count = overlap_counts.combined.all,
  annotation = "GENCODE+",
  group = "all"
)
closeAllConnections()

## Combine all data
exon_overlap_data <- rbind(
  overlap_data.assembly.all,
  overlap_data.gencode.all,
  overlap_data.combined.all
)

exon_overlap_data$annotation <- factor(exon_overlap_data$annotation,
                                       levels = c("GENCODEv45", "GENCODE+", "lr-assembly"))

## Scatter plot: InfRV vs overlapping exon count
# Merge overlap data with InfRV data
scatter_data <- exon_overlap_data %>%
  left_join(infrv.GUSTO.plot %>% dplyr::select(isoform, InfRV, annotation), 
            by = c("transcript_id" = "isoform", "annotation" = "annotation"))

# Calculate R² and p-values for each annotation group
regression_stats <- scatter_data %>%
  group_by(annotation) %>%
  do(tidy(summary(lm(log10(overlap_count+1) ~ log10(InfRV), data = .)))) %>%
  filter(term == "log10(InfRV)") %>%
  ungroup()

r_squared <- scatter_data %>%
  group_by(annotation) %>%
  do(data.frame(r2 = summary(lm(log10(overlap_count+1) ~ log10(InfRV), data = .))$r.squared)) %>%
  ungroup()

annotation_stats <- left_join(regression_stats, r_squared, by = "annotation")

r_values <- scatter_data %>%
  dplyr::group_by(annotation) %>%
  dplyr::summarise(r = cor(log10(overlap_count+1), log10(InfRV), use = "complete.obs")) %>%
  ungroup()

annotation_stats <- annotation_stats %>%
  left_join(r_values, by = "annotation")

annotation_labels <- data.frame(
  x = -2,
  y = c(3, 2.7, 2.4),
  label = sprintf("r = %.3f, R² = %.3f, p = %.2e", 
                  annotation_stats$r, 
                  annotation_stats$r2, 
                  annotation_stats$p.value),
  annotation = annotation_stats$annotation
)

figs6h <- ggplot(scatter_data, aes(x = log10(InfRV), y = log10(overlap_count+1), color = annotation)) +
  geom_point(alpha = 0.1, size = 1) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
  geom_label(data = annotation_labels, aes(x = x, y = y, label = label),
             hjust = 0, size = 3, show.legend = FALSE, 
             fill = "white", label.size = 0.5) +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank()) +
  scale_color_manual(values = c("royalblue", "#de77ae", "seagreen")) +
  xlab("log10(InfRV)") +
  ylab("log10(count of shared exons)") +
  labs(title="InfRV scales with exon sharing",
       subtitle="Data from GUSTO (n=200)")
figs6h

ggsave("Figures/figs6h.pdf", figs6h, width = 7.38, height = 4.71, units = "in", dpi = 300)

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
GUSTO.exprCorr.gencode$group <- "All"
GUSTO.exprCorr.gencode$annotation <- "GENCODEv45"

## GENCODE+ correlations
counts.GUSTO.combined <- counts(estimateSizeFactors(dds.combined.GUSTO), normalized = TRUE)
GUSTO.exprCorr.combined.FSM <- pcorr.bystruc(counts.GUSTO.combined, txi.combined.GUSTO$length, classif.pcorr, "FSM")
GUSTO.exprCorr.combined.ISM <- pcorr.bystruc(counts.GUSTO.combined, txi.combined.GUSTO$length, classif.pcorr, "ISM")
GUSTO.exprCorr.combined.NIC <- pcorr.bystruc(counts.GUSTO.combined, txi.combined.GUSTO$length, classif.pcorr, "NIC")
GUSTO.exprCorr.combined.NNC <- pcorr.bystruc(counts.GUSTO.combined, txi.combined.GUSTO$length, classif.pcorr, "NNC")
GUSTO.exprCorr.combined.Other <- pcorr.bystruc(counts.GUSTO.combined, txi.combined.GUSTO$length, classif.pcorr, "Other")

# Add annotation column to all structural categories
GUSTO.exprCorr.assembly.FSM$annotation <- "lr-assembly"
GUSTO.exprCorr.assembly.ISM$annotation <- "lr-assembly"
GUSTO.exprCorr.assembly.NIC$annotation <- "lr-assembly"
GUSTO.exprCorr.assembly.NNC$annotation <- "lr-assembly"
GUSTO.exprCorr.assembly.Other$annotation <- "lr-assembly"

GUSTO.exprCorr.combined.FSM$annotation <- "GENCODE+"
GUSTO.exprCorr.combined.ISM$annotation <- "GENCODE+"
GUSTO.exprCorr.combined.NIC$annotation <- "GENCODE+"
GUSTO.exprCorr.combined.NNC$annotation <- "GENCODE+"
GUSTO.exprCorr.combined.Other$annotation <- "GENCODE+"

# GENCODE+ all transcripts
GUSTO.datExpr.TPM.combined <- counts.GUSTO.combined / txi.combined.GUSTO$length
GUSTO.datExpr.TPM.combined <- log(t(t(GUSTO.datExpr.TPM.combined) * 1e6 / colSums(GUSTO.datExpr.TPM.combined)) + 1)
GUSTO.exprCorr.combined.all <- pcorr_parallel(GUSTO.datExpr.TPM.combined)
GUSTO.exprCorr.combined.all$group <- "All"
GUSTO.exprCorr.combined.all$annotation <- "GENCODE+"

# lr-assembly all transcripts
GUSTO.datExpr.TPM.assembly <- counts.GUSTO.assembly / txi.assembly.GUSTO$length
GUSTO.datExpr.TPM.assembly <- log(t(t(GUSTO.datExpr.TPM.assembly) * 1e6 / colSums(GUSTO.datExpr.TPM.assembly)) + 1)
GUSTO.exprCorr.assembly.all <- pcorr_parallel(GUSTO.datExpr.TPM.assembly)
GUSTO.exprCorr.assembly.all$group <- "All"
GUSTO.exprCorr.assembly.all$annotation <- "lr-assembly"

# Combine all correlation results
GUSTO.exprCorr.merged <- rbind(GUSTO.exprCorr.gencode, GUSTO.exprCorr.combined.all, GUSTO.exprCorr.assembly.all,
                               GUSTO.exprCorr.assembly.FSM, GUSTO.exprCorr.combined.FSM,
                               GUSTO.exprCorr.assembly.ISM, GUSTO.exprCorr.combined.ISM,
                               GUSTO.exprCorr.assembly.NIC, GUSTO.exprCorr.combined.NIC,
                               GUSTO.exprCorr.assembly.NNC, GUSTO.exprCorr.combined.NNC,
                               GUSTO.exprCorr.assembly.Other, GUSTO.exprCorr.combined.Other)

GUSTO.exprCorr.merged$group <- factor(GUSTO.exprCorr.merged$group,
                                      levels = c("All", "FSM", "ISM", "NIC", "NNC", "Other"))
GUSTO.exprCorr.merged$annotation <- factor(GUSTO.exprCorr.merged$annotation,
                                           levels = c("GENCODEv45", "GENCODE+", "lr-assembly"))

# Prepare data - lr-assembly only
figs6h.dat <- GUSTO.exprCorr.merged %>%
  subset(annotation == "lr-assembly")

# Create grouped categories
figs6h.dat$group_categorized <- as.character(figs6h.dat$group)
figs6h.dat$group_categorized[!figs6h.dat$group %in% c("All", "FSM", "ISM", "NIC", "NNC")] <- "Other"

# Set factor levels
figs6h.dat$group_categorized <- factor(figs6h.dat$group_categorized,
                                       levels = c("All", "FSM", "ISM", "NIC", "NNC", "Other"))


################################################################################
# FigS6F: Pairwise sample correlations in GUSTO
################################################################################

# Color palette
myPalette <- c("white", "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "black")

# Create plot
figs6f <- ggplot(figs6h.dat, aes(x = group_categorized, fill = group_categorized, y = R)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.7) +
  theme_minimal() + 
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.title.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        plot.title=element_text(size=12)) +
  scale_fill_manual(values = myPalette) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("R") + 
  ylim(0, 1) + labs(subtitle="Data from GUSTO (n=200)") +
  ggtitle("Pairwise sample correlations by structural category")

figs6f
ggsave("Figures/figs6f.pdf", figs6f, width = 4.58, height = 2.98, units = "in", dpi = 300)

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
Gen3G.exprCorr.gencode$group <- "All"
Gen3G.exprCorr.gencode$annotation <- "GENCODEv45"

## Gen3G GENCODE+ correlations
counts.Gen3G.combined <- counts(estimateSizeFactors(dds.combined.Gen3G), normalized = TRUE)
Gen3G.exprCorr.combined.FSM <- pcorr.bystruc(counts.Gen3G.combined, txi.combined.Gen3G$length, classif.pcorr, "FSM")
Gen3G.exprCorr.combined.ISM <- pcorr.bystruc(counts.Gen3G.combined, txi.combined.Gen3G$length, classif.pcorr, "ISM")
Gen3G.exprCorr.combined.NIC <- pcorr.bystruc(counts.Gen3G.combined, txi.combined.Gen3G$length, classif.pcorr, "NIC")
Gen3G.exprCorr.combined.NNC <- pcorr.bystruc(counts.Gen3G.combined, txi.combined.Gen3G$length, classif.pcorr, "NNC")
Gen3G.exprCorr.combined.Other <- pcorr.bystruc(counts.Gen3G.combined, txi.combined.Gen3G$length, classif.pcorr, "Other")

# Add annotation column to all structural categories
Gen3G.exprCorr.assembly.FSM$annotation <- "lr-assembly"
Gen3G.exprCorr.assembly.ISM$annotation <- "lr-assembly"
Gen3G.exprCorr.assembly.NIC$annotation <- "lr-assembly"
Gen3G.exprCorr.assembly.NNC$annotation <- "lr-assembly"
Gen3G.exprCorr.assembly.Other$annotation <- "lr-assembly"

Gen3G.exprCorr.combined.FSM$annotation <- "GENCODE+"
Gen3G.exprCorr.combined.ISM$annotation <- "GENCODE+"
Gen3G.exprCorr.combined.NIC$annotation <- "GENCODE+"
Gen3G.exprCorr.combined.NNC$annotation <- "GENCODE+"
Gen3G.exprCorr.combined.Other$annotation <- "GENCODE+"

# GENCODE+ all transcripts
Gen3G.datExpr.TPM.combined <- counts.Gen3G.combined / txi.combined.Gen3G$length
Gen3G.datExpr.TPM.combined <- log(t(t(Gen3G.datExpr.TPM.combined) * 1e6 / colSums(Gen3G.datExpr.TPM.combined)) + 1)
Gen3G.exprCorr.combined.all <- pcorr_parallel(Gen3G.datExpr.TPM.combined)
Gen3G.exprCorr.combined.all$group <- "All"
Gen3G.exprCorr.combined.all$annotation <- "GENCODE+"

# lr-assembly all transcripts
Gen3G.datExpr.TPM.assembly <- counts.Gen3G.assembly / txi.assembly.Gen3G$length
Gen3G.datExpr.TPM.assembly <- log(t(t(Gen3G.datExpr.TPM.assembly) * 1e6 / colSums(Gen3G.datExpr.TPM.assembly)) + 1)
Gen3G.exprCorr.assembly.all <- pcorr_parallel(Gen3G.datExpr.TPM.assembly)
Gen3G.exprCorr.assembly.all$group <- "All"
Gen3G.exprCorr.assembly.all$annotation <- "lr-assembly"

# Combine Gen3G correlation results
Gen3G.exprCorr.merged <- rbind(Gen3G.exprCorr.gencode, Gen3G.exprCorr.combined.all, Gen3G.exprCorr.assembly.all,
                               Gen3G.exprCorr.assembly.FSM, Gen3G.exprCorr.combined.FSM,
                               Gen3G.exprCorr.assembly.ISM, Gen3G.exprCorr.combined.ISM,
                               Gen3G.exprCorr.assembly.NIC, Gen3G.exprCorr.combined.NIC,
                               Gen3G.exprCorr.assembly.NNC, Gen3G.exprCorr.combined.NNC,
                               Gen3G.exprCorr.assembly.Other, Gen3G.exprCorr.combined.Other)

Gen3G.exprCorr.merged$group <- factor(Gen3G.exprCorr.merged$group,
                                      levels = c("All", "FSM", "ISM", "NIC", "NNC", "Other"))
Gen3G.exprCorr.merged$annotation <- factor(Gen3G.exprCorr.merged$annotation,
                                           levels = c("GENCODEv45", "GENCODE+", "lr-assembly"))

# Prepare data - lr-assembly only
figs6i.dat <- Gen3G.exprCorr.merged %>%
  subset(annotation == "lr-assembly")
# Create grouped categories
figs6i.dat$group_categorized <- as.character(figs6i.dat$group)
figs6i.dat$group_categorized[!figs6i.dat$group %in% c("All", "FSM", "ISM", "NIC", "NNC")] <- "Other"
# Set factor levels
figs6i.dat$group_categorized <- factor(figs6i.dat$group_categorized,
                                       levels = c("All", "FSM", "ISM", "NIC", "NNC", "Other"))

################################################################################
# FigS6G: Pairwise sample correlations in Gen3G
################################################################################

# Color palette
myPalette <- c("white", "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "black")
# Create plot
figs6g <- ggplot(figs6i.dat, aes(x = group_categorized, fill = group_categorized, y = R)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.7) +
  theme_minimal() + 
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.title.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        plot.title=element_text(size=12)) +
  scale_fill_manual(values = myPalette) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("R") + 
  ylim(0, 1) + labs(subtitle="Data from Gen3G (n=152)") +
  ggtitle("Pairwise sample correlations by structural category")
figs6g

ggsave("Figures/figs6g.pdf", figs6i, width = 4.58, height = 2.98, units = "in", dpi = 300)

################################################################################
# End of Analysis
################################################################################