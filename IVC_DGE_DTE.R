################################################################################
# Differential Gene and Transcript Expression Analysis: GDM Association Study
# Author: Sean T. Bresnahan
# Description: Comprehensive differential expression analysis comparing gestational
#              diabetes mellitus (GDM) cases vs controls using both gene-level
#              (DEG) and transcript-level (DTE) approaches. Analyzes expression
#              patterns using filtered placental transcriptome assembly versus
#              GENCODE v45 reference across two independent cohorts.
#
# Input Datasets:
#   - Assembly/ESPRESSO_corrected_SC_filtered_classification.txt
#       SQANTI3 structural classification for filtered assembly
#   - Supplementary/tx2g_ESPRESSO_assembly_SC_filtered.RData
#       Transcript-to-gene mapping objects (tx2g.assembly, tx2g.gencode, tx2g.combined)
#   - Supplementary/metadata_GUSTO.RData
#       GUSTO cohort metadata (metadata.GUSTO, cordat.GUSTO)
#   - Supplementary/metadata_Gen3G.RData
#       Gen3G cohort metadata (metadata.Gen3G, cordat.Gen3G.Gen3G)
#   - /rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/infrv_GUSTO/COUNTS_assembly_SC/*/quant.sf
#       Salmon quantification files for GUSTO cohort using lr-assembly
#   - /rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/infrv_GUSTO/COUNTS_gencode/*/quant.sf
#       Salmon quantification files for GUSTO cohort using GENCODE v45
#   - /rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/infrv_GUSTO/COUNTS_combined/*/quant.sf
#       Salmon quantification files for GUSTO cohort using GENCODE+ (combined)
#   - /rsrch5/home/epi/stbresnahan/scratch/natcomms_revision/COUNTS_GUSTO_Adipose_Subcutaneous/*/quant.sf
#       Salmon quantification for GUSTO samples using GTEx Adipose assembly
#   - /rsrch5/home/epi/stbresnahan/scratch/natcomms_revision/COUNTS_GUSTO_Cells_Cultured_fibroblasts/*/quant.sf
#       Salmon quantification for GUSTO samples using GTEx Fibroblasts assembly
#   - /rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/infrv_Gen3G/COUNTS_assembly_SC/*/quant.sf
#       Salmon quantification files for Gen3G cohort using lr-assembly
#   - /rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/infrv_Gen3G/COUNTS_gencode/*/quant.sf
#       Salmon quantification files for Gen3G cohort using GENCODE v45
#   - /rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/infrv_Gen3G/COUNTS_combined/*/quant.sf
#       Salmon quantification files for Gen3G cohort using GENCODE+ (combined)
#   - /rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/combined_GTF/lr_assembly_GENCODE_v45_combined.gtf
#       Combined GENCODE v45 + filtered assembly GTF (GENCODE+)
#   - /rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/GTEx_v9/SQANTI3/Adipose_Subcutaneous_filtered/Adipose_Subcutaneous_filtered.fixed.gtf
#       GTEx Adipose tissue assembly GTF
#   - /rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/GTEx_v9/SQANTI3/Cells_Cultured_fibroblasts_filtered/Cells_Cultured_fibroblasts_filtered.fixed.gtf
#       GTEx Fibroblasts assembly GTF
#   - DIR_ALIGN/*_gencode.tsv.gz
#       Read-to-transcript mapping files for GENCODE v45 (per sample)
#   - DIR_ALIGN/*_assembly.tsv.gz
#       Read-to-transcript mapping files for lr-assembly (per sample)
#
# Output Files:
#   Data Files:
#     - Results/GUSTO_ESPRESSO_assembly_filtered_SC.RData
#         GUSTO differential expression results (all annotations and tissues)
#     - Results/Gen3G_ESPRESSO_assembly_filtered_SC.RData
#         Gen3G differential expression results (all annotations)
#     - Supplementary/TableS10.csv
#         Differential transcript expression results across cohorts and annotations
#     - Supplementary/TableS11.csv
#         Differential gene expression results across cohorts and annotations
#
#   Figures - Differential Expression Overview:
#     - Figures/figs7a.pdf
#         Volcano plots for DTE across cohorts and annotations (GUSTO, Gen3G)
#     - Figures/figs7b.pdf
#         Volcano plots for DEG across cohorts and annotations (GUSTO, Gen3G)
#     - Figures/fig3e.pdf
#         UpSet plot comparing DTE and DEG across cohorts and annotations
#     - Figures/figs8c.pdf
#         UpSet plot for GUSTO tissues (Placenta, Adipose, Fibroblasts)
#
#   Figures - CSH1 Gene Analysis:
#     - Figures/fig5b.pdf
#         Log fold change comparison for CSH1 FSM transcripts (lr-assembly vs GENCODEv45)
#     - Figures/figs5c.pdf
#         Read count correlation for CSH1 (ENST00000558284.1 highlighted)
#     - Figures/fig5d.pdf
#         Alluvial diagram showing CSH1 read mapping transitions (GENCODEv45 → lr-assembly)
#
#   Figures - GNAS Gene Analysis:
#     - Figures/figs5e.pdf
#         Log fold change comparison for GNAS FSM transcripts (lr-assembly vs GENCODEv45)
#     - Figures/figs5f.pdf
#         Read count correlation for GNAS (ENST00000476196.5 highlighted)
#     - Figures/fig5g.pdf
#         Alluvial diagram showing GNAS read mapping transitions (GENCODEv45 → lr-assembly)
#
################################################################################

################################################################################
# Environment Setup
################################################################################

# Configure R library paths for HPC environment
.libPaths(c("/home/stbresnahan/R/ubuntu/4.3.1", .libPaths()))

# Load required libraries
library(dplyr)         # Data manipulation
library(rtracklayer)   # Genomic data import/export
library(ggplot2)       # Data visualization
library(edgeR)         # RNA-seq analysis toolkit
library(tidyr)         # Data tidying
library(DESeq2)        # Differential expression analysis
library(emmeans)       # Estimated marginal means
library(deming)        # Deming regression
library(GenomicFeatures) # Genomic feature manipulation
library(pcaExplorer)   # PCA analysis and visualization
library(plyr)          # Data manipulation (legacy functions)
library(PupillometryR) # Statistical functions
library(WGCNA)         # Weighted gene co-expression network analysis
library(pheatmap)      # Heatmap visualization
library(org.Hs.eg.db)  # Human genome annotation
library(clusterProfiler) # Gene enrichment analysis
library(doParallel)    # Parallel processing backend
library(doSNOW)        # Snow parallel backend
library(ggfortify)     # Statistical visualization extensions
library(tximport)      # Transcript abundance import and summarization
library(sva)           # Surrogate variable analysis
library(CompDTUReg)    # Comparative differential transcript usage
library(RUVSeq)        # Remove unwanted variation
library(SummarizedExperiment) # Bioconductor data containers
library(isotwas)       # Isoform-level TWAS methods
library(BiocParallel)  # Bioconductor parallel processing
library(ggExtra)       # Additional ggplot2 functionality
library(fishpond)      # Transcript-level inferential variance
library(tximeta)       # Metadata-aware transcript quantification import
library(cowplot)       # Plot arrangement
library(scales)        # Scale functions for ggplot2
library(ggdist)        # Distribution visualization
library(ggpubr)        # Publication-ready plots
library(tibble)        # Enhanced data frames
library(ggVennDiagram) # Venn diagram visualization
library(ggalluvial)    # Alluvial diagram visualization
library(data.table)    # High-performance data manipulation
library(readxl)        # Excel file reading
library(ggupset)

################################################################################
# Global Configuration
################################################################################

# Define color palettes for consistent visualization
myPalette <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#ffd92f",
               "#e5c494", "#c87774", "#d3adaf", "#b3b3b3")
myPalette2 <- c(myPalette[1:4], "grey20")

################################################################################
# Core Functions
################################################################################

# DESeq2 analysis with RUVSeq batch correction
deseq2_ruv <- function(SumEx, meta, tx = FALSE, tx2g, trait,
                       trait.type = c("categorical", "continuous"),
                       type.cat.levels, sample.ID.column, batch, weights = TRUE, parallel = TRUE, W = 1) {
  return.list <- list()
  
  # Standardize column names
  names(meta)[which(names(meta) == trait)] <- "trait"
  if(exists("batch")) {names(meta)[which(names(meta) == batch)] <- "batch"}
  
  # Filter low-expressed genes and remove NA samples
  counts <- SumEx$counts[rowSums(SumEx$counts >= 5) >= 3, ]
  counts <- counts[, !colnames(counts) %in% meta[is.na(meta$trait), sample.ID.column]]
  meta <- meta[!is.na(meta$trait), ]
  
  ### RUVr batch correction
  print("Making DGEList")
  design <- model.matrix(~trait, data = meta)
  y <- DGEList(counts = counts, group = meta$trait)
  print("calcNormFactors")
  y <- calcNormFactors(y, method = "upperquartile", na.rm = TRUE)
  print("estimateGLMCommonDisp")
  y <- estimateGLMCommonDisp(y, design)
  print("estimateGLMTagwiseDisp")
  y <- estimateGLMTagwiseDisp(y, design)
  
  print("Fitting model and returning residuals")
  fit <- glmFit(y, design)
  res <- counts - fit$fitted.values
  print("Performing RUVSeq")
  set4 <- RUVSeq::RUVr(x = round(counts), residuals = res, 
                       k = 2, cIdx = rownames(counts))
  
  # Optional batch effect visualization
  if(exists("batch")) {
    pca <- prcomp(t(counts))
    pca.pc1 <- data.frame(PC1 = pca[["x"]][, 1],
                          Batch = meta$batch)
    boxplot(PC1 ~ Batch, pca.pc1)
  }
  
  pca.ruv <- prcomp(t(set4[["normalizedCounts"]]))
  
  if(exists("batch")) {
    pca.ruv.pc1 <- data.frame(PC1 = pca.ruv[["x"]][, 1],
                              Batch = meta$batch)
    boxplot(PC1 ~ Batch, pca.ruv.pc1)
  }
  
  return.list[[1]] <- pca.ruv
  
  # Add RUV factors to metadata
  meta$W1 <- set4$W[, 1]
  meta$W2 <- set4$W[, 2]
  
  return.list[[2]] <- meta
  
  ### DESeq2 analysis with RUV adjustment
  register(MulticoreParam(4))
  
  print("Performing DESeq2")
  if(weights == TRUE) {
    if(W == 1) {
      dds <- DESeqDataSetFromTximport(txi = SumEx, colData = meta, 
                                      design = ~ sex + GA + W1 + trait)
    } else {
      dds <- DESeqDataSetFromTximport(txi = SumEx, colData = meta, 
                                      design = ~ sex + GA + W1 + W2 + trait)
    }
  } else {
    dds <- DESeqDataSetFromTximport(txi = SumEx, colData = meta, 
                                    design = ~ sex + GA + trait)
  }
  
  dds <- DESeq(dds, parallel = parallel)
  
  # Extract results based on trait type
  if(trait.type == "categorical") {
    res <- data.frame(results(dds, contrast = c("trait", type.cat.levels)))
  }
  if(trait.type == "continuous") {
    res <- data.frame(results(dds, name = "trait"))
  }
  
  # Clean results
  res <- res[!is.na(res$pvalue), ]
  res <- res[!is.na(res$padj), ]
  
  hist(res$pvalue)
  
  return.list[[3]] <- dds
  return.list[[4]] <- res
  
  # Optional transcript-level gene aggregation
  if(tx) {
    print("tx = TRUE, running gene-level aggregation")
    res$transcript_id <- row.names(res)
    res <- left_join(res, tx2g[, c("transcript_id", "gene_id")], by = "transcript_id")
    
    # Screen for significant genes
    gene_res <- res %>%
      group_by(gene_id) %>%
      summarise(Screen.P = isotwas::p_screen(pvalue))
    gene_res$Screen.P.Adj <- p.adjust(gene_res$Screen.P, 'BH')
    
    ttt <- left_join(res, gene_res, by = 'gene_id')
    
    # Calculate confirmation threshold
    alpha2 <- (sum(gene_res$Screen.P.Adj < 0.05) * 0.05) / length(unique(ttt$gene_id))
    
    # Confirm significant isoforms within significant genes
    isoform_res <- ttt %>%
      group_by(gene_id) %>%
      summarise(transcript_id = transcript_id,
                baseMean = baseMean,
                log2FoldChange = log2FoldChange,
                lfcSE = lfcSE,
                stat = stat,
                pvalue = pvalue,
                Screen.P.Adj = Screen.P.Adj,
                Confirmation.P = isotwas::p_confirm(pvalue, alpha = alpha2))
    
    return.list[[5]] <- isoform_res
  }
  
  return.list[[6]] <- set4
  
  return(return.list)
}

# Function to filter tximport object by sample IDs
filter_tximport <- function(x, filter.list) {
  x1 <- data.frame(x[["abundance"]])
  x[["abundance"]] <- as.matrix(x1[, !names(x1) %in% filter.list])
  x2 <- data.frame(x[["counts"]])
  x[["counts"]] <- as.matrix(x2[, !names(x2) %in% filter.list])
  x3 <- data.frame(x[["length"]])
  x[["length"]] <- as.matrix(x3[, !names(x3) %in% filter.list])
  return(x)
}

################################################################################
# Data Import and Processing
################################################################################

# Assembly isoform classification from SQANTI3
classif <- read.table("Assembly/ESPRESSO_corrected_SC_filtered_classification.txt", 
                      header = TRUE)
names(classif)[1] <- "transcript_id"

# Standardize structural category definitions
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

################################################################################
# Transcript-to-Gene Mapping tables
################################################################################

load("Supplementary/tx2g_ESPRESSO_assembly_SC_filtered.RData")

### Calculate isoforms per gene for assembly
iso.tab.assembly <- data.frame(table(tx2g.assembly$gene_id))
names(iso.tab.assembly) <- c("gene_id", "Freq.assembly")
iso.tab.assembly$gene_id <- as.character(iso.tab.assembly$gene_id)

## Calculate isoforms per gene for GENCODE
iso.tab.gencode <- data.frame(table(tx2g.gencode$gene_id))
names(iso.tab.gencode) <- c("gene_id", "Freq.gencode")
iso.tab.gencode$gene_id <- as.character(iso.tab.gencode$gene_id)

################################################################################
# GUSTO Cohort Analysis
################################################################################

## Import GUSTO metadata
load("Supplementary/metadata_GUSTO.RData")

################################################################################
# GUSTO Assembly Analysis
################################################################################

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

# Define samples to filter (outliers + missing GDM data)
n <- colnames(txi.assembly.g[["counts"]])
filter.assembly <- c(n[!n%in%cordat.GUSTO$SampleID])
txi.assembly.filter <- filter_tximport(txi.assembly.g, filter.assembly)

### Gene-level differential expression analysis
drgC.assembly.gdm.gene <- deseq2_ruv(txi.assembly.filter, cordat.GUSTO, tx = FALSE, tx2g.assembly,
                                     "gdm", "categorical", c("0", "1"), "SampleID",
                                     "Batch", weights = TRUE, W = 2)

## Transcript-level expression analysis with QU correction
s.assembly <- catchSalmon(dirs.salmon.assembly)

### Import transcript-level quantification with QU correction
txi.assembly <- tximport(unlist(files.list.assembly), type = "salmon", tx2gene = tx2g.assembly,
                         geneIdCol = "gene_id", txIdCol = "transcript_id", dropInfReps = TRUE, txOut = TRUE)
txi.assembly$counts <- txi.assembly$counts / s.assembly$annotation$Overdispersion
txi.assembly$abundance <- cpm(txi.assembly$counts)

### Filter and analyze transcript-level expression
txi.assembly.filter <- filter_tximport(txi.assembly, filter.assembly)

#### Differential transcript expression for GDM
drgC.assembly.gdm <- deseq2_ruv(txi.assembly.filter, cordat.GUSTO, tx = FALSE, tx2g.assembly,
                                "gdm", "categorical", c("0", "1"), "SampleID",
                                "Batch", weights = TRUE, W = 2)

################################################################################
# GUSTO GENCODE Analysis
################################################################################

## Gene-level expression analysis with GENCODE
dirs.salmon.gencode <- list.dirs("/rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/infrv_GUSTO/COUNTS_gencode", recursive = FALSE)
files.list.gencode <- list()
for (i in 1:length(dirs.salmon.gencode)) {
  files.list.gencode[i] <- paste(unlist(dirs.salmon.gencode[i]), "quant.sf", sep = "/")
}
dirnames.gencode <- sapply(strsplit(unlist(dirs.salmon.gencode), "/"), "[",
                           length(unlist(strsplit(unlist(dirs.salmon.gencode), "/")[1])))
names(files.list.gencode) <- dirnames.gencode

# Import GENCODE gene-level quantification
txi.gencode.g <- tximport(unlist(files.list.gencode), type = "salmon", tx2gene = tx2g.gencode,
                          geneIdCol = "gene_id", txIdCol = "transcript_id", dropInfReps = TRUE, txOut = FALSE)

# Apply same filtering as assembly analysis
filter.gencode <- filter.assembly
txi.gencode.filter <- filter_tximport(txi.gencode.g, filter.gencode)

# Gene-level differential expression
drgC.gencode.gdm.gene <- deseq2_ruv(txi.gencode.filter, cordat.GUSTO, tx = FALSE, tx2g.gencode,
                                    "gdm", "categorical", c("0", "1"), "SampleID",
                                    "Batch", weights = TRUE, W = 2)

## Transcript-level expression analysis with GENCODE
s.gencode <- catchSalmon(dirs.salmon.gencode)

### QU correction for GENCODE
txi.gencode <- tximport(unlist(files.list.gencode), type = "salmon", tx2gene = tx2g.gencode,
                        geneIdCol = "gene_id", txIdCol = "transcript_id", dropInfReps = TRUE, txOut = TRUE)
txi.gencode$counts <- txi.gencode$counts / s.gencode$annotation$Overdispersion
txi.gencode$abundance <- cpm(txi.gencode$counts)

# Filter and analyze
txi.gencode.filter <- filter_tximport(txi.gencode, filter.gencode)

### Differential transcript expression with GENCODE
drgC.gencode.gdm <- deseq2_ruv(txi.gencode.filter, cordat.GUSTO, tx = FALSE, tx2g.gencode,
                               "gdm", "categorical", c("0", "1"), "SampleID",
                               "Batch", weights = TRUE, W = 2)

################################################################################
# GUSTO GENCODE+ Analysis
################################################################################
tx2g.combined <- data.frame(import("/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/combined_GTF/lr_assembly_GENCODE_v45_combined.gtf"))
tx2g.combined <- tx2g.combined[tx2g.combined$type=="transcript",c("transcript_id","gene_id")]
tx2g.combined <- tx2g.combined[!duplicated(tx2g.combined),]

## Gene-level expression analysis with GENCODE_PLUS
dirs.salmon.gencode_plus <- list.dirs("/rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/infrv_GUSTO/COUNTS_combined", recursive = FALSE)
files.list.gencode_plus <- list()
for (i in 1:length(dirs.salmon.gencode_plus)) {
  files.list.gencode_plus[i] <- paste(unlist(dirs.salmon.gencode_plus[i]), "quant.sf", sep = "/")
}
dirnames.gencode_plus <- sapply(strsplit(unlist(dirs.salmon.gencode_plus), "/"), "[",
                                length(unlist(strsplit(unlist(dirs.salmon.gencode_plus), "/")[1])))
names(files.list.gencode_plus) <- dirnames.gencode_plus
# Import GENCODE_PLUS gene-level quantification
txi.gencode_plus.g <- tximport(unlist(files.list.gencode_plus), type = "salmon", tx2gene = tx2g.combined,
                               geneIdCol = "gene_id", txIdCol = "transcript_id", dropInfReps = TRUE, txOut = FALSE)
# Apply same filtering as assembly analysis
filter.gencode_plus <- filter.assembly
txi.gencode_plus.filter <- filter_tximport(txi.gencode_plus.g, filter.gencode_plus)
# Gene-level differential expression
drgC.gencode_plus.gdm.gene <- deseq2_ruv(txi.gencode_plus.filter, cordat.GUSTO, tx = FALSE, tx2g.combined,
                                         "gdm", "categorical", c("0", "1"), "SampleID",
                                         "Batch", weights = TRUE, W = 2)
## Transcript-level expression analysis with GENCODE_PLUS
s.gencode_plus <- catchSalmon(dirs.salmon.gencode_plus)
### QU correction for GENCODE_PLUS
txi.gencode_plus <- tximport(unlist(files.list.gencode_plus), type = "salmon", tx2gene = tx2g.combined,
                             geneIdCol = "gene_id", txIdCol = "transcript_id", dropInfReps = TRUE, txOut = TRUE)
txi.gencode_plus$counts <- txi.gencode_plus$counts / s.gencode_plus$annotation$Overdispersion
txi.gencode_plus$abundance <- cpm(txi.gencode_plus$counts)
# Filter and analyze
txi.gencode_plus.filter <- filter_tximport(txi.gencode_plus, filter.gencode_plus)
### Differential transcript expression with GENCODE_PLUS
drgC.gencode_plus.gdm <- deseq2_ruv(txi.gencode_plus.filter, cordat.GUSTO, tx = FALSE, tx2g.combined,
                                    "gdm", "categorical", c("0", "1"), "SampleID",
                                    "Batch", weights = TRUE, W = 2)

################################################################################
# GUSTO Adipose Analysis
################################################################################
tx2g.adipose <- data.frame(import("/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/GTEx_v9/SQANTI3/Adipose_Subcutaneous_filtered/Adipose_Subcutaneous_filtered.fixed.gtf"))
tx2g.adipose <- tx2g.adipose[tx2g.adipose$type=="transcript",c("transcript_id","gene_id")]
tx2g.adipose <- tx2g.adipose[!duplicated(tx2g.adipose),]
## Gene-level expression analysis with Adipose
dirs.salmon.adipose <- list.dirs("/rsrch5/home/epi/stbresnahan/scratch/natcomms_revision/COUNTS_GUSTO_Adipose_Subcutaneous", recursive = FALSE)
files.list.adipose <- list()
for (i in 1:length(dirs.salmon.adipose)) {
  files.list.adipose[i] <- paste(unlist(dirs.salmon.adipose[i]), "quant.sf", sep = "/")
}
dirnames.adipose <- sapply(strsplit(unlist(dirs.salmon.adipose), "/"), "[",
                           length(unlist(strsplit(unlist(dirs.salmon.adipose), "/")[1])))
names(files.list.adipose) <- dirnames.adipose
# Import Adipose gene-level quantification
txi.adipose.g <- tximport(unlist(files.list.adipose), type = "salmon", tx2gene = tx2g.adipose,
                          geneIdCol = "gene_id", txIdCol = "transcript_id", dropInfReps = TRUE, txOut = FALSE)
# Apply same filtering as assembly analysis
filter.adipose <- filter.assembly
txi.adipose.filter <- filter_tximport(txi.adipose.g, filter.adipose)
# Gene-level differential expression
drgC.adipose.gdm.gene <- deseq2_ruv(txi.adipose.filter, cordat.GUSTO, tx = FALSE, tx2g.adipose,
                                    "gdm", "categorical", c("0", "1"), "SampleID",
                                    "Batch", weights = TRUE, W = 2)
## Transcript-level expression analysis with Adipose
s.adipose <- catchSalmon(dirs.salmon.adipose)
### QU correction for Adipose
txi.adipose <- tximport(unlist(files.list.adipose), type = "salmon", tx2gene = tx2g.adipose,
                        geneIdCol = "gene_id", txIdCol = "transcript_id", dropInfReps = TRUE, txOut = TRUE)
txi.adipose$counts <- txi.adipose$counts / s.adipose$annotation$Overdispersion
txi.adipose$abundance <- cpm(txi.adipose$counts)
# Filter and analyze
txi.adipose.filter <- filter_tximport(txi.adipose, filter.adipose)
### Differential transcript expression with Adipose
drgC.adipose.gdm <- deseq2_ruv(txi.adipose.filter, cordat.GUSTO, tx = FALSE, tx2g.adipose,
                               "gdm", "categorical", c("0", "1"), "SampleID",
                               "Batch", weights = TRUE, W = 2)

################################################################################
# GUSTO Fibroblasts Analysis
################################################################################
tx2g.fibroblasts <- data.frame(import("/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/GTEx_v9/SQANTI3/Cells_Cultured_fibroblasts_filtered/Cells_Cultured_fibroblasts_filtered.fixed.gtf"))
tx2g.fibroblasts <- tx2g.fibroblasts[tx2g.fibroblasts$type=="transcript",c("transcript_id","gene_id")]
tx2g.fibroblasts <- tx2g.fibroblasts[!duplicated(tx2g.fibroblasts),]
## Gene-level expression analysis with Fibroblasts
dirs.salmon.fibroblasts <- list.dirs("/rsrch5/home/epi/stbresnahan/scratch/natcomms_revision/COUNTS_GUSTO_Cells_Cultured_fibroblasts", recursive = FALSE)
files.list.fibroblasts <- list()
for (i in 1:length(dirs.salmon.fibroblasts)) {
  files.list.fibroblasts[i] <- paste(unlist(dirs.salmon.fibroblasts[i]), "quant.sf", sep = "/")
}
dirnames.fibroblasts <- sapply(strsplit(unlist(dirs.salmon.fibroblasts), "/"), "[",
                               length(unlist(strsplit(unlist(dirs.salmon.fibroblasts), "/")[1])))
names(files.list.fibroblasts) <- dirnames.fibroblasts
# Import Fibroblasts gene-level quantification
txi.fibroblasts.g <- tximport(unlist(files.list.fibroblasts), type = "salmon", tx2gene = tx2g.fibroblasts,
                              geneIdCol = "gene_id", txIdCol = "transcript_id", dropInfReps = TRUE, txOut = FALSE)
# Apply same filtering as assembly analysis
filter.fibroblasts <- filter.assembly
txi.fibroblasts.filter <- filter_tximport(txi.fibroblasts.g, filter.fibroblasts)
# Gene-level differential expression
drgC.fibroblasts.gdm.gene <- deseq2_ruv(txi.fibroblasts.filter, cordat.GUSTO, tx = FALSE, tx2g.fibroblasts,
                                        "gdm", "categorical", c("0", "1"), "SampleID",
                                        "Batch", weights = TRUE, W = 2)
## Transcript-level expression analysis with Fibroblasts
s.fibroblasts <- catchSalmon(dirs.salmon.fibroblasts)
### QU correction for Fibroblasts
txi.fibroblasts <- tximport(unlist(files.list.fibroblasts), type = "salmon", tx2gene = tx2g.fibroblasts,
                            geneIdCol = "gene_id", txIdCol = "transcript_id", dropInfReps = TRUE, txOut = TRUE)
txi.fibroblasts$counts <- txi.fibroblasts$counts / s.fibroblasts$annotation$Overdispersion
txi.fibroblasts$abundance <- cpm(txi.fibroblasts$counts)
# Filter and analyze
txi.fibroblasts.filter <- filter_tximport(txi.fibroblasts, filter.fibroblasts)
### Differential transcript expression with Fibroblasts
drgC.fibroblasts.gdm <- deseq2_ruv(txi.fibroblasts.filter, cordat.GUSTO, tx = FALSE, tx2g.fibroblasts,
                                   "gdm", "categorical", c("0", "1"), "SampleID",
                                   "Batch", weights = TRUE, W = 2)

################################################################################
# Save GUSTO results for cross-cohort comparison
################################################################################
GUSTO.drgC.assembly.gdm <- drgC.assembly.gdm
GUSTO.drgC.assembly.gdm.gene <- drgC.assembly.gdm.gene
GUSTO.drgC.gencode.gdm <- drgC.gencode.gdm
GUSTO.drgC.gencode.gdm.gene <- drgC.gencode.gdm.gene
GUSTO.drgC.gencode_plus.gdm <- drgC.gencode_plus.gdm
GUSTO.drgC.gencode_plus.gdm.gene <- drgC.gencode_plus.gdm.gene
GUSTO.drgC.adipose.gdm <- drgC.adipose.gdm
GUSTO.drgC.adipose.gdm.gene <- drgC.adipose.gdm.gene
GUSTO.drgC.fibroblasts.gdm <- drgC.fibroblasts.gdm
GUSTO.drgC.fibroblasts.gdm.gene <- drgC.fibroblasts.gdm.gene
save(list = c("GUSTO.drgC.assembly.gdm",
              "GUSTO.drgC.gencode.gdm",
              "GUSTO.drgC.gencode_plus.gdm",
              "GUSTO.drgC.assembly.gdm.gene",
              "GUSTO.drgC.gencode.gdm.gene",
              "GUSTO.drgC.gencode_plus.gdm.gene",
              "GUSTO.drgC.adipose.gdm",
              "GUSTO.drgC.adipose.gdm.gene",
              "GUSTO.drgC.fibroblasts.gdm",
              "GUSTO.drgC.fibroblasts.gdm.gene"),
     file = "Results/GUSTO_ESPRESSO_assembly_filtered_SC.RData")

################################################################################
# Gen3G Cohort Analysis
################################################################################

## Import Gen3G metadata
load("Supplementary/metadata_Gen3G.RData")

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

# Filter out problematic samples and bad sequencing batches
n <- colnames(txi.assembly.g[["counts"]])
filter.assembly <- c(n[!n%in%cordat.Gen3G$Run])
txi.assembly.filter <- filter_tximport(txi.assembly.g, filter.assembly)

# Gene-level differential expression
drgC.assembly.gdm.gene <- deseq2_ruv(txi.assembly.filter, cordat.Gen3G, tx = FALSE, tx2g.assembly,
                                     "gdm", "categorical", c("0", "1"), "Run",
                                     "sequencing_batch", weights = TRUE, W = 2)

## Transcript-level expression analysis with QU correction
s.assembly <- catchSalmon(dirs.salmon.assembly)

### QU correction for assembly
txi.assembly <- tximport(unlist(files.list.assembly), type = "salmon", tx2gene = tx2g.assembly,
                         geneIdCol = "gene_id", txIdCol = "transcript_id", dropInfReps = TRUE, txOut = TRUE)
txi.assembly$counts <- txi.assembly$counts / s.assembly$annotation$Overdispersion
txi.assembly$abundance <- cpm(txi.assembly$counts)

### Filter and analyze transcript-level expression
txi.assembly.filter <- filter_tximport(txi.assembly, filter.assembly)

#### Differential transcript expression for GDM
drgC.assembly.gdm <- deseq2_ruv(txi.assembly.filter, cordat.Gen3G, tx = F, tx2g.assembly,
                                "gdm", "categorical", c("0", "1"), "Run",
                                "sequencing_batch", weights = TRUE, W = 2)

################################################################################
# Gen3G GENCODE Analysis
################################################################################

## Prepare GENCODE gene expression data
dirs.salmon.gencode <- list.dirs("/rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/infrv_Gen3G/COUNTS_gencode", 
                                 recursive = FALSE)
files.list.gencode <- list()
for (i in 1:length(dirs.salmon.gencode)) {
  files.list.gencode[i] <- paste(unlist(dirs.salmon.gencode[i]), "quant.sf", sep = "/")
}
dirnames.gencode <- sapply(strsplit(unlist(dirs.salmon.gencode), "/"), "[",
                           length(unlist(strsplit(unlist(dirs.salmon.gencode), "/")[1])))
names(files.list.gencode) <- dirnames.gencode

# Filter to samples present in metadata
keep <- which(dirnames.gencode %in% cordat.Gen3G$Run)
dirs.salmon.gencode <- dirs.salmon.gencode[keep]
dirnames.gencode <- dirnames.gencode[keep]
files.list.gencode <- unlist(files.list.gencode)[keep]

# Import GENCODE gene-level quantification
txi.gencode.g <- tximport(unlist(files.list.gencode), type = "salmon", tx2gene = tx2g.gencode,
                          geneIdCol = "gene_id", txIdCol = "transcript_id", dropInfReps = TRUE, txOut = FALSE)

# Apply same filtering as assembly analysis
filter.gencode <- filter.assembly
txi.gencode.filter <- filter_tximport(txi.gencode.g, filter.gencode)

# Gene-level differential expression
drgC.gencode.gdm.gene <- deseq2_ruv(txi.gencode.filter, cordat.Gen3G, tx = FALSE, tx2g.gencode,
                                    "gdm", "categorical", c("0", "1"), "Run",
                                    "sequencing_batch", weights = TRUE, W = 2)

## Transcript-level expression analysis with GENCODE
s.gencode <- catchSalmon(dirs.salmon.gencode)

### QU correction for GENCODE
txi.gencode <- tximport(unlist(files.list.gencode), type = "salmon", tx2gene = tx2g.gencode,
                        geneIdCol = "gene_id", txIdCol = "transcript_id", dropInfReps = TRUE, txOut = TRUE)
txi.gencode$counts <- txi.gencode$counts / s.gencode$annotation$Overdispersion
txi.gencode$abundance <- cpm(txi.gencode$counts)

# Filter and analyze
txi.gencode.filter <- filter_tximport(txi.gencode, filter.gencode)

#### Differential transcript expression with GENCODE
drgC.gencode.gdm <- deseq2_ruv(txi.gencode.filter, cordat.Gen3G, tx = F, tx2g.gencode,
                               "gdm", "categorical", c("0", "1"), "Run",
                               "sequencing_batch", weights = TRUE, W = 2)

################################################################################
# Gen3G GENCODE+ Analysis
################################################################################
## Prepare GENCODE_PLUS gene expression data
dirs.salmon.gencode_plus <- list.dirs("/rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/infrv_Gen3G/COUNTS_combined", 
                                      recursive = FALSE)
files.list.gencode_plus <- list()
for (i in 1:length(dirs.salmon.gencode_plus)) {
  files.list.gencode_plus[i] <- paste(unlist(dirs.salmon.gencode_plus[i]), "quant.sf", sep = "/")
}
dirnames.gencode_plus <- sapply(strsplit(unlist(dirs.salmon.gencode_plus), "/"), "[",
                                length(unlist(strsplit(unlist(dirs.salmon.gencode_plus), "/")[1])))
names(files.list.gencode_plus) <- dirnames.gencode_plus
# Filter to samples present in metadata
keep <- which(dirnames.gencode_plus %in% cordat.Gen3G$Run)
dirs.salmon.gencode_plus <- dirs.salmon.gencode_plus[keep]
dirnames.gencode_plus <- dirnames.gencode_plus[keep]
files.list.gencode_plus <- unlist(files.list.gencode_plus)[keep]
# Import GENCODE_PLUS gene-level quantification
txi.gencode_plus.g <- tximport(unlist(files.list.gencode_plus), type = "salmon", tx2gene = tx2g.combined,
                               geneIdCol = "gene_id", txIdCol = "transcript_id", dropInfReps = TRUE, txOut = FALSE)
# Apply same filtering as assembly analysis
filter.gencode_plus <- filter.assembly
txi.gencode_plus.filter <- filter_tximport(txi.gencode_plus.g, filter.gencode_plus)
# Gene-level differential expression
drgC.gencode_plus.gdm.gene <- deseq2_ruv(txi.gencode_plus.filter, cordat.Gen3G, tx = FALSE, tx2g.combined,
                                         "gdm", "categorical", c("0", "1"), "Run",
                                         "sequencing_batch", weights = TRUE, W = 2)
## Transcript-level expression analysis with GENCODE_PLUS
s.gencode_plus <- catchSalmon(dirs.salmon.gencode_plus)
### QU correction for GENCODE_PLUS
txi.gencode_plus <- tximport(unlist(files.list.gencode_plus), type = "salmon", tx2gene = tx2g.combined,
                             geneIdCol = "gene_id", txIdCol = "transcript_id", dropInfReps = TRUE, txOut = TRUE)
txi.gencode_plus$counts <- txi.gencode_plus$counts / s.gencode_plus$annotation$Overdispersion
txi.gencode_plus$abundance <- cpm(txi.gencode_plus$counts)
# Filter and analyze
txi.gencode_plus.filter <- filter_tximport(txi.gencode_plus, filter.gencode_plus)
#### Differential transcript expression with GENCODE_PLUS
drgC.gencode_plus.gdm <- deseq2_ruv(txi.gencode_plus.filter, cordat.Gen3G, tx = F, tx2g.combined,
                                    "gdm", "categorical", c("0", "1"), "Run",
                                    "sequencing_batch", weights = TRUE, W = 2)

################################################################################
# Save Gen3G results for cross-cohort comparison
################################################################################
Gen3G.drgC.assembly.gdm <- drgC.assembly.gdm
Gen3G.drgC.assembly.gdm.gene <- drgC.assembly.gdm.gene
Gen3G.drgC.gencode.gdm <- drgC.gencode.gdm
Gen3G.drgC.gencode.gdm.gene <- drgC.gencode.gdm.gene
Gen3G.drgC.gencode_plus.gdm <- drgC.gencode_plus.gdm
Gen3G.drgC.gencode_plus.gdm.gene <- drgC.gencode_plus.gdm.gene
save(list = c("Gen3G.drgC.assembly.gdm", "Gen3G.drgC.gencode.gdm", "Gen3G.drgC.gencode_plus.gdm",
              "Gen3G.drgC.assembly.gdm.gene", "Gen3G.drgC.gencode.gdm.gene", "Gen3G.drgC.gencode_plus.gdm.gene"),
     file = "Results/Gen3G_ESPRESSO_assembly_filtered_SC.RData")




################################################################################
# Cross-Cohort Comparative Analysis
################################################################################

# Prepare classification data for visualization
classif.pcorr <- classif
classif.pcorr$structural_category <- revalue(classif.pcorr$structural_category, 
                                             c("Genic Genomic" = "Other",
                                               "Antisense" = "Other",
                                               "Fusion" = "Other",
                                               "Intergenic" = "Other",
                                               "Genic Intron" = "Other"))

################################################################################
# VST Transformed Data for Visualization
################################################################################

## VST transformed read counts for visualization
GUSTO.vst.assembly <- assay(varianceStabilizingTransformation(GUSTO.drgC.assembly.gdm[[3]], blind = FALSE))
GUSTO.vst.gencode <- assay(varianceStabilizingTransformation(GUSTO.drgC.gencode.gdm[[3]], blind = FALSE))
GUSTO.vst.gencode_plus <- assay(varianceStabilizingTransformation(GUSTO.drgC.gencode_plus.gdm[[3]], blind = FALSE))
Gen3G.vst.assembly <- assay(varianceStabilizingTransformation(Gen3G.drgC.assembly.gdm[[3]], blind = FALSE))
Gen3G.vst.gencode <- assay(varianceStabilizingTransformation(Gen3G.drgC.gencode.gdm[[3]], blind = FALSE))
Gen3G.vst.gencode_plus <- assay(varianceStabilizingTransformation(Gen3G.drgC.gencode_plus.gdm[[3]], blind = FALSE))

################################################################################
# Results Tables Preparation
################################################################################

## Results tables for differential expression
results.assembly.GUSTO.gdm <- GUSTO.drgC.assembly.gdm[[4]]
results.assembly.GUSTO.gdm$transcript_id <- row.names(results.assembly.GUSTO.gdm)
results.assembly.Gen3G.gdm <- Gen3G.drgC.assembly.gdm[[4]]
results.assembly.Gen3G.gdm$transcript_id <- row.names(results.assembly.Gen3G.gdm)

results.gencode.GUSTO.gdm <- GUSTO.drgC.gencode.gdm[[4]]
results.gencode.GUSTO.gdm$transcript_id <- row.names(results.gencode.GUSTO.gdm)
results.gencode.Gen3G.gdm <- Gen3G.drgC.gencode.gdm[[4]]
results.gencode.Gen3G.gdm$transcript_id <- row.names(results.gencode.Gen3G.gdm)

results.gencode_plus.GUSTO.gdm <- GUSTO.drgC.gencode_plus.gdm[[4]]
results.gencode_plus.GUSTO.gdm$transcript_id <- row.names(results.gencode_plus.GUSTO.gdm)
results.gencode_plus.Gen3G.gdm <- Gen3G.drgC.gencode_plus.gdm[[4]]
results.gencode_plus.Gen3G.gdm$transcript_id <- row.names(results.gencode_plus.Gen3G.gdm)

results.assembly.GUSTO.gdm.g <- GUSTO.drgC.assembly.gdm[[5]]
results.assembly.Gen3G.gdm.g <- Gen3G.drgC.assembly.gdm[[5]]
results.gencode.GUSTO.gdm.g <- GUSTO.drgC.gencode.gdm[[5]]
results.gencode.Gen3G.gdm.g <- Gen3G.drgC.gencode.gdm[[5]]
results.gencode_plus.GUSTO.gdm.g <- GUSTO.drgC.gencode_plus.gdm[[5]]
results.gencode_plus.Gen3G.gdm.g <- Gen3G.drgC.gencode_plus.gdm[[5]]

results.assembly.GUSTO.gdm.gene <- GUSTO.drgC.assembly.gdm.gene[[4]]
results.assembly.GUSTO.gdm.gene$gene_id <- row.names(results.assembly.GUSTO.gdm.gene)
results.assembly.Gen3G.gdm.gene <- Gen3G.drgC.assembly.gdm.gene[[4]]
results.assembly.Gen3G.gdm.gene$gene_id <- row.names(results.assembly.Gen3G.gdm.gene)

results.gencode.GUSTO.gdm.gene <- GUSTO.drgC.gencode.gdm.gene[[4]]
results.gencode.GUSTO.gdm.gene$gene_id <- row.names(results.gencode.GUSTO.gdm.gene)
results.gencode.Gen3G.gdm.gene <- Gen3G.drgC.gencode.gdm.gene[[4]]
results.gencode.Gen3G.gdm.gene$gene_id <- row.names(results.gencode.Gen3G.gdm.gene)

results.gencode_plus.GUSTO.gdm.gene <- GUSTO.drgC.gencode_plus.gdm.gene[[4]]
results.gencode_plus.GUSTO.gdm.gene$gene_id <- row.names(results.gencode_plus.GUSTO.gdm.gene)
results.gencode_plus.Gen3G.gdm.gene <- Gen3G.drgC.gencode_plus.gdm.gene[[4]]
results.gencode_plus.Gen3G.gdm.gene$gene_id <- row.names(results.gencode_plus.Gen3G.gdm.gene)

################################################################################
# Significant Gene/Transcript Lists
################################################################################

## Total numbers of DTEs and DGEs at FDR < 0.1
GUSTO.gdm.DTE.assembly <- results.assembly.GUSTO.gdm[results.assembly.GUSTO.gdm$padj < 0.1, "transcript_id"]
GUSTO.gdm.DTE.gencode <- results.gencode.GUSTO.gdm[results.gencode.GUSTO.gdm$padj < 0.1, "transcript_id"]
GUSTO.gdm.DTE.gencode_plus <- results.gencode_plus.GUSTO.gdm[results.gencode_plus.GUSTO.gdm$padj < 0.1, "transcript_id"]
Gen3G.gdm.DTE.assembly <- results.assembly.Gen3G.gdm[results.assembly.Gen3G.gdm$padj < 0.1, "transcript_id"]
Gen3G.gdm.DTE.gencode <- results.gencode.Gen3G.gdm[results.gencode.Gen3G.gdm$padj < 0.1, "transcript_id"]
Gen3G.gdm.DTE.gencode_plus <- results.gencode_plus.Gen3G.gdm[results.gencode_plus.Gen3G.gdm$padj < 0.1, "transcript_id"]

GUSTO.gdm.DGE.assembly <- results.assembly.GUSTO.gdm.gene[results.assembly.GUSTO.gdm.gene$padj < 0.1, "gene_id"]
GUSTO.gdm.DGE.gencode <- results.gencode.GUSTO.gdm.gene[results.gencode.GUSTO.gdm.gene$padj < 0.1, "gene_id"]
GUSTO.gdm.DGE.gencode_plus <- results.gencode_plus.GUSTO.gdm.gene[results.gencode_plus.GUSTO.gdm.gene$padj < 0.1, "gene_id"]
Gen3G.gdm.DGE.assembly <- results.assembly.Gen3G.gdm.gene[results.assembly.Gen3G.gdm.gene$padj < 0.1, "gene_id"]
Gen3G.gdm.DGE.gencode <- results.gencode.Gen3G.gdm.gene[results.gencode.Gen3G.gdm.gene$padj < 0.1, "gene_id"]
Gen3G.gdm.DGE.gencode_plus <- results.gencode_plus.Gen3G.gdm.gene[results.gencode_plus.Gen3G.gdm.gene$padj < 0.1, "gene_id"]

################################################################################
# Function to thin points in large plots
################################################################################

thin_points <- function(df, x_col = "log2FoldChange", y_col = "neglog10_padj",
                        sig_col = "sig", bins = 200, max_points_per_bin = 5) {
  
  library(dplyr)
  
  # Separate significant and non-significant points
  sig_points <- df %>% filter(!!sym(sig_col))
  non_sig <- df %>% filter(! (!!sym(sig_col)))
  
  # Compute bin coordinates
  x_breaks <- seq(min(non_sig[[x_col]]), max(non_sig[[x_col]]), length.out = bins + 1)
  y_breaks <- seq(min(non_sig[[y_col]]), max(non_sig[[y_col]]), length.out = bins + 1)
  
  # Assign each point to a bin
  non_sig <- non_sig %>%
    mutate(
      x_bin = cut(.data[[x_col]], breaks = x_breaks, include.lowest = TRUE, labels = FALSE),
      y_bin = cut(.data[[y_col]], breaks = y_breaks, include.lowest = TRUE, labels = FALSE)
    )
  
  # Sample points per bin safely
  non_sig_thinned <- non_sig %>%
    group_by(x_bin, y_bin) %>%
    slice_head(n = max_points_per_bin) %>%  # take first N points per bin
    ungroup() %>%
    select(-x_bin, -y_bin)
  
  # Combine with significant points
  bind_rows(sig_points, non_sig_thinned)
}

################################################################################
# Figure S7A: DET volcano plots
################################################################################
results.assembly.GUSTO.gdm$Cohort <- "GUSTO"
results.assembly.GUSTO.gdm$Annotation <- "lr-assembly"

results.gencode.GUSTO.gdm$Cohort <- "GUSTO"
results.gencode.GUSTO.gdm$Annotation <- "GENCODEv45"

results.gencode_plus.GUSTO.gdm$Cohort <- "GUSTO"
results.gencode_plus.GUSTO.gdm$Annotation <- "GENCODE+"

results.assembly.Gen3G.gdm$Cohort <- "Gen3G"
results.assembly.Gen3G.gdm$Annotation <- "lr-assembly"

results.gencode.Gen3G.gdm$Cohort <- "Gen3G"
results.gencode.Gen3G.gdm$Annotation <- "GENCODEv45"

results.gencode_plus.Gen3G.gdm$Cohort <- "Gen3G"
results.gencode_plus.Gen3G.gdm$Annotation <- "GENCODE+"

dat.fig7a <- rbind(results.assembly.GUSTO.gdm,
                   results.gencode.GUSTO.gdm,
                   results.gencode_plus.GUSTO.gdm,
                   results.assembly.Gen3G.gdm,
                   results.gencode.Gen3G.gdm,
                   results.gencode_plus.Gen3G.gdm)

write.csv(dat.fig7a,"Supplementary/TableS10.csv")

dat.fig7a <- dat.fig7a %>%
  mutate(
    neglog10_padj = -log10(padj),
    sig = padj < 0.1
  )

dat.fig7a.thinned <- thin_points(dat.fig7a,
                           x_col = "log2FoldChange",
                           y_col = "neglog10_padj",
                           sig_col = "sig",
                           bins = 200,
                           max_points_per_bin = 5)


figs7a <- ggplot(dat.fig7a.thinned,
                    aes(x = log2FoldChange, y = neglog10_padj, color = sig)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 12) +
  labs(
    x = "log2 Fold Change",
    y = "-log10(adjusted p-value)",
    color = "Significant"
  ) +
  facet_grid(Cohort ~ Annotation) +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    panel.grid.minor = element_blank()
  ) + ggtitle("Differential transcript expression by GDM status")

figs7a

ggsave("Figures/figs7a.pdf", figs7a, width = 7.12, height = 5.5, units = "in", dpi = 300)

################################################################################
# Figure S7B: DEG volcano plots
################################################################################
results.assembly.GUSTO.gdm.gene$Cohort <- "GUSTO"
results.assembly.GUSTO.gdm.gene$Annotation <- "lr-assembly"

results.gencode.GUSTO.gdm.gene$Cohort <- "GUSTO"
results.gencode.GUSTO.gdm.gene$Annotation <- "GENCODEv45"

results.gencode_plus.GUSTO.gdm.gene$Cohort <- "GUSTO"
results.gencode_plus.GUSTO.gdm.gene$Annotation <- "GENCODE+"

results.assembly.Gen3G.gdm.gene$Cohort <- "Gen3G"
results.assembly.Gen3G.gdm.gene$Annotation <- "lr-assembly"

results.gencode.Gen3G.gdm.gene$Cohort <- "Gen3G"
results.gencode.Gen3G.gdm.gene$Annotation <- "GENCODEv45"

results.gencode_plus.Gen3G.gdm.gene$Cohort <- "Gen3G"
results.gencode_plus.Gen3G.gdm.gene$Annotation <- "GENCODE+"

dat.fig7b <- rbind(results.assembly.GUSTO.gdm.gene,
                   results.gencode.GUSTO.gdm.gene,
                   results.gencode_plus.GUSTO.gdm.gene,
                   results.assembly.Gen3G.gdm.gene,
                   results.gencode.Gen3G.gdm.gene,
                   results.gencode_plus.Gen3G.gdm.gene)

write.csv(dat.fig7b,"Supplementary/TableS11.csv")

dat.fig7b <- dat.fig7b %>%
  mutate(
    neglog10_padj = -log10(padj),
    sig = padj < 0.1
  )

dat.fig7b.thinned <- thin_points(dat.fig7b,
                                 x_col = "log2FoldChange",
                                 y_col = "neglog10_padj",
                                 sig_col = "sig",
                                 bins = 200,
                                 max_points_per_bin = 5)


figs7b <- ggplot(dat.fig7b.thinned,
                    aes(x = log2FoldChange, y = neglog10_padj, color = sig)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 12) +
  labs(
    x = "log2 Fold Change",
    y = "-log10(adjusted p-value)",
    color = "Significant"
  ) +
  facet_grid(Cohort ~ Annotation) +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    panel.grid.minor = element_blank()
  ) + ggtitle("Differential gene expression by GDM status")

ggsave("Figures/figs7b.pdf", figs7b, width = 7.12, height = 5.5, units = "in", dpi = 300)


################################################################################
# Figure 3E: UpSet Plot Comparing DTE and DGE
################################################################################

# -----------------------------
# Prepare data for UpSet plot
# -----------------------------
upset_data <- data.frame(
  gene = unique(c(
    GUSTO.gdm.DTE.assembly,
    GUSTO.gdm.DTE.gencode,
    GUSTO.gdm.DTE.gencode_plus,
    Gen3G.gdm.DTE.assembly,
    Gen3G.gdm.DTE.gencode,
    Gen3G.gdm.DTE.gencode_plus,
    GUSTO.gdm.DGE.assembly,
    GUSTO.gdm.DGE.gencode,
    GUSTO.gdm.DGE.gencode_plus,
    Gen3G.gdm.DGE.assembly,
    Gen3G.gdm.DGE.gencode,
    Gen3G.gdm.DGE.gencode_plus
  ))
)

upset_data$sets <- lapply(upset_data$gene, function(g) {
  sets <- c()
  if (g %in% GUSTO.gdm.DTE.assembly) sets <- c(sets, "GUSTO lr-assembly DET")
  if (g %in% GUSTO.gdm.DTE.gencode) sets <- c(sets, "GUSTO GENCODEv45 DET")
  if (g %in% GUSTO.gdm.DTE.gencode_plus) sets <- c(sets, "GUSTO GENCODE+ DET")
  if (g %in% Gen3G.gdm.DTE.assembly) sets <- c(sets, "Gen3G lr-assembly DET")
  if (g %in% Gen3G.gdm.DTE.gencode) sets <- c(sets, "Gen3G GENCODEv45 DET")
  if (g %in% Gen3G.gdm.DTE.gencode_plus) sets <- c(sets, "Gen3G GENCODE+ DET")
  if (g %in% GUSTO.gdm.DGE.assembly) sets <- c(sets, "GUSTO lr-assembly DEG")
  if (g %in% GUSTO.gdm.DGE.gencode) sets <- c(sets, "GUSTO GENCODEv45 DEG")
  if (g %in% GUSTO.gdm.DGE.gencode_plus) sets <- c(sets, "GUSTO GENCODE+ DEG")
  if (g %in% Gen3G.gdm.DGE.assembly) sets <- c(sets, "Gen3G lr-assembly DEG")
  if (g %in% Gen3G.gdm.DGE.gencode) sets <- c(sets, "Gen3G GENCODEv45 DEG")
  if (g %in% Gen3G.gdm.DGE.gencode_plus) sets <- c(sets, "Gen3G GENCODE+ DEG")
  return(sets)
})

# -----------------------------
# Filter to intersections >= 10
# -----------------------------
upset_counts <- table(sapply(upset_data$sets, function(x) paste(sort(x), collapse = "|")))
upset_filtered <- upset_data[sapply(upset_data$sets, function(x) {
  key <- paste(sort(x), collapse = "|")
  upset_counts[key] >= 10
}), ]

# -----------------------------
# Compute set order and colors
# -----------------------------
set_order <- upset_filtered %>%
  tidyr::unnest_longer(sets) %>%
  dplyr::count(sets, name = "freq") %>%
  dplyr::arrange(desc(freq), sets) %>%
  dplyr::pull(sets)

# Assign colors by pattern (do NOT reverse here)
set_colors <- sapply(set_order, function(set_name) {
  if (grepl("lr-assembly", set_name)) return("seagreen")
  if (grepl("GENCODEv45", set_name)) return("royalblue")
  if (grepl("GENCODE\\+", set_name)) return("#de77ae")
  return("black")
})

set_colors_ordered <- unname(set_colors)  # for matching later

# -----------------------------
# Plot UpSet figure
# -----------------------------
fig3e <- ggplot(upset_filtered, aes(x = sets)) +
  geom_bar() +
  scale_x_upset(n_intersections = Inf) +
  theme_minimal() +
  ylab("Count") +
  labs(
    title = "Differential expression across cohorts and annotations",
    subtitle = expression("Displaying sets of n" >= "10 features")
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = -100)), 
    plot.margin = margin(l = 100, t = 10),
    plot.title.position = "plot",       # title relative to entire plot area
    plot.title = element_text(hjust = 0),     
    plot.subtitle = element_text(hjust = 0)
  ) +
  axis_combmatrix(
    override_plotting_function = function(df) {
      # Force top-to-bottom order
      df$single_label <- factor(df$single_label, levels = rev(set_order))
      
      # Assign point and label colors
      df <- df %>%
        dplyr::mutate(
          point_color = ifelse(
            observed,
            set_colors_ordered[match(as.character(single_label), set_order)],
            "grey95"
          ),
          label_color = set_colors_ordered[match(as.character(single_label), set_order)]
        )
      
      # Plot combination matrix
      ggplot(df, aes(x = at, y = single_label)) +
        geom_point(aes(color = point_color), size = 3) +
        geom_line(
          data = function(dat) dat[dat$observed, , drop = FALSE],
          aes(group = labels),
          color = "grey30",
          linewidth = 1.2
        ) +
        # Colored labels
        scale_y_discrete(
          labels = function(labs) {
            lab_cols <- set_colors_ordered[match(labs, set_order)]
            purrr::map2_chr(labs, lab_cols, ~ paste0("<span style='color:", .y, "'>", .x, "</span>"))
          }
        ) +
        ylab("") + xlab("") +
        scale_color_identity() +
        scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
        theme_void() +
        theme(
          axis.text.y = ggtext::element_markdown(size = 10, hjust = 1),
          axis.text.x = element_blank()
        )
    }
  )

fig3e

ggsave("Figures/fig3e.pdf", fig3e, width = 5.76, height = 3.68, units = "in", dpi = 300)



################################################################################
# Additional Results Tables for Adipose and Fibroblasts
################################################################################

results.adipose.GUSTO.gdm <- GUSTO.drgC.adipose.gdm[[4]]
results.adipose.GUSTO.gdm$transcript_id <- row.names(results.adipose.GUSTO.gdm)
results.fibroblasts.GUSTO.gdm <- GUSTO.drgC.fibroblasts.gdm[[4]]
results.fibroblasts.GUSTO.gdm$transcript_id <- row.names(results.fibroblasts.GUSTO.gdm)

results.adipose.GUSTO.gdm.gene <- GUSTO.drgC.adipose.gdm.gene[[4]]
results.adipose.GUSTO.gdm.gene$gene_id <- row.names(results.adipose.GUSTO.gdm.gene)
results.fibroblasts.GUSTO.gdm.gene <- GUSTO.drgC.fibroblasts.gdm.gene[[4]]
results.fibroblasts.GUSTO.gdm.gene$gene_id <- row.names(results.fibroblasts.GUSTO.gdm.gene)

################################################################################
# Significant Gene/Transcript Lists for Adipose and Fibroblasts
################################################################################

GUSTO.gdm.DTE.adipose <- results.adipose.GUSTO.gdm[results.adipose.GUSTO.gdm$padj < 0.1, "transcript_id"]
GUSTO.gdm.DTE.fibroblasts <- results.fibroblasts.GUSTO.gdm[results.fibroblasts.GUSTO.gdm$padj < 0.1, "transcript_id"]
GUSTO.gdm.DGE.adipose <- results.adipose.GUSTO.gdm.gene[results.adipose.GUSTO.gdm.gene$padj < 0.1, "gene_id"]
GUSTO.gdm.DGE.fibroblasts <- results.fibroblasts.GUSTO.gdm.gene[results.fibroblasts.GUSTO.gdm.gene$padj < 0.1, "gene_id"]

################################################################################
# Figure S8C: UpSet Plot for GUSTO Tissues (Including Adipose and Fibroblasts)
################################################################################

# -----------------------------
# Prepare data for GUSTO-only UpSet plot
# -----------------------------
upset_data_gusto <- data.frame(
  gene = unique(c(
    GUSTO.gdm.DTE.assembly,
    GUSTO.gdm.DTE.gencode,
    GUSTO.gdm.DTE.gencode_plus,
    GUSTO.gdm.DTE.adipose,
    GUSTO.gdm.DTE.fibroblasts,
    GUSTO.gdm.DGE.assembly,
    GUSTO.gdm.DGE.gencode,
    GUSTO.gdm.DGE.gencode_plus,
    GUSTO.gdm.DGE.adipose,
    GUSTO.gdm.DGE.fibroblasts
  ))
)

upset_data_gusto$sets <- lapply(upset_data_gusto$gene, function(g) {
  sets <- c()
  if (g %in% GUSTO.gdm.DTE.assembly) sets <- c(sets, "Placenta lr-assembly DET")
  if (g %in% GUSTO.gdm.DTE.gencode) sets <- c(sets, "Placenta GENCODEv45 DET")
  if (g %in% GUSTO.gdm.DTE.gencode_plus) sets <- c(sets, "Placenta GENCODE+ DET")
  if (g %in% GUSTO.gdm.DTE.adipose) sets <- c(sets, "Adipose DET")
  if (g %in% GUSTO.gdm.DTE.fibroblasts) sets <- c(sets, "Fibroblasts DET")
  if (g %in% GUSTO.gdm.DGE.assembly) sets <- c(sets, "Placenta lr-assembly DEG")
  if (g %in% GUSTO.gdm.DGE.gencode) sets <- c(sets, "Placenta GENCODEv45 DEG")
  if (g %in% GUSTO.gdm.DGE.gencode_plus) sets <- c(sets, "Placenta GENCODE+ DEG")
  if (g %in% GUSTO.gdm.DGE.adipose) sets <- c(sets, "Adipose DEG")
  if (g %in% GUSTO.gdm.DGE.fibroblasts) sets <- c(sets, "Fibroblasts DEG")
  return(sets)
})

# -----------------------------
# Filter to intersections >= 10
# -----------------------------
upset_counts_gusto <- table(sapply(upset_data_gusto$sets, function(x) paste(sort(x), collapse = "|")))
upset_filtered_gusto <- upset_data_gusto[sapply(upset_data_gusto$sets, function(x) {
  key <- paste(sort(x), collapse = "|")
  upset_counts_gusto[key] >= 10
}), ]

# -----------------------------
# Compute set order and colors
# -----------------------------
set_order_gusto <- upset_filtered_gusto %>%
  tidyr::unnest_longer(sets) %>%
  dplyr::count(sets, name = "freq") %>%
  dplyr::arrange(desc(freq), sets) %>%
  dplyr::pull(sets)

# Assign colors by pattern
set_colors_gusto <- sapply(set_order_gusto, function(set_name) {
  if (grepl("Adipose", set_name)) return("#ff68a1")
  if (grepl("Fibroblasts", set_name)) return("#e68613")
  if (grepl("lr-assembly", set_name)) return("seagreen")
  if (grepl("GENCODEv45", set_name)) return("royalblue")
  if (grepl("GENCODE\\+", set_name)) return("#de77ae")
  return("black")
})

set_colors_ordered_gusto <- unname(set_colors_gusto)

# -----------------------------
# Plot UpSet figure for GUSTO tissues
# -----------------------------

figs8c <- ggplot(upset_filtered_gusto, aes(x = sets)) +
  geom_bar() +
  scale_x_upset(n_intersections = Inf) +
  theme_minimal() +
  ylab("Count") +
  labs(
    title = "Differential expression in GUSTO across annotations",
    subtitle = expression("Displaying sets of n" >= "10 features")
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = -100)), 
    plot.margin = margin(l = 100, t = 10),
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0),     
    plot.subtitle = element_text(hjust = 0)
  ) +
  axis_combmatrix(
    override_plotting_function = function(df) {
      # Force top-to-bottom order
      df$single_label <- factor(df$single_label, levels = rev(set_order_gusto))
      
      # Assign point and label colors
      df <- df %>%
        dplyr::mutate(
          point_color = ifelse(
            observed,
            set_colors_ordered_gusto[match(as.character(single_label), set_order_gusto)],
            "grey95"
          ),
          label_color = set_colors_ordered_gusto[match(as.character(single_label), set_order_gusto)]
        )
      
      # Plot combination matrix
      ggplot(df, aes(x = at, y = single_label)) +
        geom_point(aes(color = point_color), size = 3) +
        geom_line(
          data = function(dat) dat[dat$observed, , drop = FALSE],
          aes(group = labels),
          color = "grey30",
          linewidth = 1.2
        ) +
        # Colored labels
        scale_y_discrete(
          labels = function(labs) {
            lab_cols <- set_colors_ordered_gusto[match(labs, set_order_gusto)]
            purrr::map2_chr(labs, lab_cols, ~ paste0("<span style='color:", .y, "'>", .x, "</span>"))
          }
        ) +
        ylab("") + xlab("") +
        scale_color_identity() +
        scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
        theme_void() +
        theme(
          axis.text.y = ggtext::element_markdown(size = 10, hjust = 1),
          axis.text.x = element_blank()
        )
    }
  )

figs8c

ggsave("Figures/figs8c.pdf", figs8c, width = 5.76, height = 3.68, units = "in", dpi = 300)


################################################################################
# Figure 5B: CSH1 Log Fold Change Comparison
################################################################################

## CSH1 full-splice match transcripts
CSH1 <- c("ENST00000453363.7", "ENST00000316193.13",
          "ENST00000329882.8", "ENST00000558284.1",
          "ENST00000610991.1", "ENST00000558661.1")

# Extract log fold changes and standard errors for assembly
dat.lfc.csh1.a <- results.assembly.GUSTO.gdm[, c("transcript_id", "log2FoldChange", "lfcSE", "pvalue")]
dat.lfc.csh1.a <- dat.lfc.csh1.a[dat.lfc.csh1.a$transcript_id %in% CSH1, ]
names(dat.lfc.csh1.a)[2:4] <- c("lfc.assembly", "lfcSE.assembly", "pvalue.assembly")

# Extract log fold changes and standard errors for GENCODE
dat.lfc.csh1.g <- results.gencode.GUSTO.gdm[, c("transcript_id", "log2FoldChange", "lfcSE", "pvalue")]
dat.lfc.csh1.g <- dat.lfc.csh1.g[dat.lfc.csh1.g$transcript_id %in% CSH1, ]
names(dat.lfc.csh1.g)[2:4] <- c("lfc.gencode", "lfcSE.gencode", "pvalue.gencode")

# Combine data and highlight specific transcript
dat.lfc.csh1 <- left_join(dat.lfc.csh1.a, dat.lfc.csh1.g)
dat.lfc.csh1 <- dat.lfc.csh1[, c(1, 5:6, 2:3)]
dat.lfc.csh1$col <- "A"
dat.lfc.csh1[4, "col"] <- "B"

highlighted <- dat.lfc.csh1[dat.lfc.csh1$col == "B", ]
non_highlighted <- dat.lfc.csh1[dat.lfc.csh1$col == "A", ]

## Figure 3B: Log fold change comparison for CSH1
fig5b <- ggplot(dat.lfc.csh1, aes(x = lfc.gencode, y = lfc.assembly)) +
  theme_classic() +
  xlab("logFC (GENCODEv45)") +
  ylab("logFC (lr-assembly)") +
  geom_errorbar(
    data = highlighted,
    aes(x = lfc.gencode,
        ymin = lfc.assembly - (lfcSE.assembly * 1.96),
        ymax = lfc.assembly + (lfcSE.assembly * 1.96)), 
    width = 0, color = "seagreen") +
  geom_errorbar(
    data = highlighted,
    aes(y = lfc.assembly,
        xmin = lfc.gencode - (lfcSE.gencode * 1.96),
        xmax = lfc.gencode + (lfcSE.gencode * 1.96)),
    width = 0, color = "royalblue") +
  geom_errorbar(
    data = non_highlighted,
    aes(x = lfc.gencode,
        ymin = lfc.assembly - (lfcSE.assembly * 1.96),
        ymax = lfc.assembly + (lfcSE.assembly * 1.96)), 
    width = 0, color = "grey80") +
  geom_errorbar(
    data = non_highlighted,
    aes(y = lfc.assembly,
        xmin = lfc.gencode - (lfcSE.gencode * 1.96),
        xmax = lfc.gencode + (lfcSE.gencode * 1.96)),
    width = 0, color = "grey70") +
  ggtitle("DTE of CSH1 FSM") +
  theme(panel.grid.major.x = element_line(color = "grey95"),
        panel.grid.major.y = element_line(color = "grey95"),
        legend.position = "none") +
  geom_hline(yintercept = 0, linetype = 'dashed', col = 'black') +
  geom_vline(xintercept = 0, linetype = 'dashed', col = 'black') +
  geom_point(data = highlighted, color = "black") +
  geom_point(data = non_highlighted, color = "grey50") +
  geom_label(data = highlighted,
             aes(label = c("ENST00000558284.1")),
             color = "black",
             fill = "white",
             label.size = 0.3,
             size = 3,
             label.r = unit(0.1, "lines"),
             nudge_y = 0.5, nudge_x = 0.3) 

fig5b
ggsave("Figures/fig5b.pdf", fig5b, width = 3, height = 3, units = "in", dpi = 300)

################################################################################
# Figure 5C: CSH1 Read Count Correlation
################################################################################

## Extract normalized read counts for CSH1 transcripts
counts.assembly.GUSTO <- data.frame(counts(GUSTO.drgC.assembly.gdm[[3]]))
counts.assembly.GUSTO <- counts.assembly.GUSTO[row.names(counts.assembly.GUSTO) %in% CSH1, ]
counts.assembly.GUSTO <- counts.assembly.GUSTO %>%
  rownames_to_column(var = "row_name") %>%
  pivot_longer(
    cols = -row_name,
    names_to = "column_name",
    values_to = "value")
names(counts.assembly.GUSTO) <- c("transcript_id", "SampleID", "counts.assembly")

counts.gencode.GUSTO <- data.frame(counts(GUSTO.drgC.gencode.gdm[[3]]))
counts.gencode.GUSTO <- counts.gencode.GUSTO[row.names(counts.gencode.GUSTO) %in% CSH1, ]
counts.gencode.GUSTO <- counts.gencode.GUSTO %>%
  rownames_to_column(var = "row_name") %>%
  pivot_longer(
    cols = -row_name,
    names_to = "column_name",
    values_to = "value")
names(counts.gencode.GUSTO) <- c("transcript_id", "SampleID", "counts.gencode")

# Combine count data and add GDM status
counts.GUSTO.CSH1 <- counts.assembly.GUSTO
counts.GUSTO.CSH1$counts.gencode <- counts.gencode.GUSTO$counts.gencode
counts.GUSTO.CSH1 <- data.frame(counts.GUSTO.CSH1)

counts.GUSTO.CSH1$col <- "A"
counts.GUSTO.CSH1[counts.GUSTO.CSH1$transcript_id == "ENST00000558284.1", "col"] <- "B"

counts.GUSTO.CSH1$col <- factor(counts.GUSTO.CSH1$col, levels = c("A", "B"))
counts.GUSTO.CSH1 <- counts.GUSTO.CSH1[order(counts.GUSTO.CSH1$col), ]
counts.GUSTO.CSH1 <- left_join(counts.GUSTO.CSH1, cordat.GUSTO[, c("SampleID", "gdm")])

counts.GUSTO.CSH1$col <- as.character(counts.GUSTO.CSH1$col)
counts.GUSTO.CSH1$gdm <- as.character(counts.GUSTO.CSH1$gdm)
counts.GUSTO.CSH1[counts.GUSTO.CSH1$gdm == 0 & counts.GUSTO.CSH1$col == "B", "col"] <- "C"

counts.GUSTO.CSH1$col <- factor(counts.GUSTO.CSH1$col, levels = c("A", "C", "B"))
counts.GUSTO.CSH1 <- counts.GUSTO.CSH1[order(counts.GUSTO.CSH1$col), ]

## Figure 3C: Read count correlation for CSH1
fig5c <- ggplot(counts.GUSTO.CSH1, aes(x = log2(counts.gencode + 1),
                                       y = log2(counts.assembly + 1),
                                       color = col)) +
  geom_point() + 
  theme_classic() +
  scale_color_manual(values = c("blue", "red", "grey50"),
                     limits = c("B", "C"),
                     name = "GDM",
                     labels = c("1", "0")) +
  xlab("GENCODEv45 counts (log2)") +
  ylab("lr-assembly counts (log2)") +
  ggtitle("Reads mapped to CSH1 FSM\n(ENST00000558284.1 highlighted)") +
  theme(panel.grid.major.x = element_line(color = "grey95"),
        panel.grid.major.y = element_line(color = "grey95"),
        legend.position = "right") +
  geom_abline(slope = 1,
              intercept = 0,
              color = "black",
              linetype = "dashed") +
  xlim(0, 25) + ylim(0, 25)

fig5c
ggsave("Figures/figs5c.pdf", fig5c, width = 3.2, height = 3, units = "in", dpi = 300)

################################################################################
# Figure 5E: GNAS Log Fold Change Comparison
################################################################################

## GNAS full-splice match transcripts
GNAS <- tx2g.gencode[tx2g.gencode$gene_id == "ENSG00000087460.29", "transcript_id"]
GNAS <- GNAS[GNAS %in% tx2g.assembly$transcript_id]

# Extract log fold changes for GNAS
dat.lfc.GNAS.a <- results.assembly.GUSTO.gdm[, c("transcript_id", "log2FoldChange", "lfcSE", "pvalue")]
dat.lfc.GNAS.a <- dat.lfc.GNAS.a[dat.lfc.GNAS.a$transcript_id %in% GNAS, ]
names(dat.lfc.GNAS.a)[2:4] <- c("lfc.assembly", "lfcSE.assembly", "pvalue.assembly")

dat.lfc.GNAS.g <- results.gencode.GUSTO.gdm[, c("transcript_id", "log2FoldChange", "lfcSE", "pvalue")]
dat.lfc.GNAS.g <- dat.lfc.GNAS.g[dat.lfc.GNAS.g$transcript_id %in% GNAS, ]
names(dat.lfc.GNAS.g)[2:4] <- c("lfc.gencode", "lfcSE.gencode", "pvalue.gencode")

# Combine and highlight specific transcript
dat.lfc.GNAS <- left_join(dat.lfc.GNAS.a, dat.lfc.GNAS.g)
dat.lfc.GNAS <- dat.lfc.GNAS[, c(1, 5:6, 2:3)]
dat.lfc.GNAS$col <- "A"
dat.lfc.GNAS[16, "col"] <- "B"

highlighted <- dat.lfc.GNAS[dat.lfc.GNAS$col == "B", ]
non_highlighted <- dat.lfc.GNAS[dat.lfc.GNAS$col == "A", ]

## Figure 3E: Log fold change comparison for GNAS
fig5e <- ggplot(dat.lfc.GNAS, aes(x = lfc.gencode, y = lfc.assembly)) +
  theme_classic() +
  xlab("logFC (GENCODEv45)") +
  ylab("logFC (lr-assembly)") +
  geom_errorbar(
    data = non_highlighted,
    aes(x = lfc.gencode,
        ymin = lfc.assembly - (lfcSE.assembly * 1.96),
        ymax = lfc.assembly + (lfcSE.assembly * 1.96)), 
    width = 0, color = "grey80") +
  geom_errorbar(
    data = non_highlighted,
    aes(y = lfc.assembly,
        xmin = lfc.gencode - (lfcSE.gencode * 1.96),
        xmax = lfc.gencode + (lfcSE.gencode * 1.96)),
    width = 0, color = "grey70") +
  geom_errorbar(
    data = highlighted,
    aes(x = lfc.gencode,
        ymin = lfc.assembly - (lfcSE.assembly * 1.96),
        ymax = lfc.assembly + (lfcSE.assembly * 1.96)), 
    width = 0, color = "seagreen") +
  geom_errorbar(
    data = highlighted,
    aes(y = lfc.assembly,
        xmin = lfc.gencode - (lfcSE.gencode * 1.96),
        xmax = lfc.gencode + (lfcSE.gencode * 1.96)),
    width = 0, color = "royalblue") +
  theme(panel.grid.major.x = element_line(color = "grey95"),
        panel.grid.major.y = element_line(color = "grey95"),
        legend.position = "none") +
  geom_hline(yintercept = 0, linetype = 'dashed', col = 'black') +
  geom_vline(xintercept = 0, linetype = 'dashed', col = 'black') +
  geom_point(data = highlighted, color = "black") +
  geom_point(data = non_highlighted, color = "grey50") +
  ggtitle("DTE of GNAS FSM") +
  geom_label(data = highlighted,
             aes(label = c("ENST00000476196.5")),
             color = "black",
             fill = "white",
             label.size = 0.3,
             size = 3,
             label.r = unit(0.1, "lines"),
             nudge_y = 0.35, nudge_x = 0) 

fig5e
ggsave("Figures/figs5e.pdf", fig5e, width = 3, height = 3, units = "in", dpi = 300)

################################################################################
# Figure 5F: GNAS Read Count Correlation
################################################################################

## Extract normalized read counts for GNAS transcripts
counts.assembly.GUSTO <- data.frame(counts(GUSTO.drgC.assembly.gdm[[3]]))
counts.assembly.GUSTO <- counts.assembly.GUSTO[row.names(counts.assembly.GUSTO) %in% GNAS, ]
counts.assembly.GUSTO <- counts.assembly.GUSTO %>%
  rownames_to_column(var = "row_name") %>%
  pivot_longer(
    cols = -row_name,
    names_to = "column_name",
    values_to = "value")
names(counts.assembly.GUSTO) <- c("transcript_id", "SampleID", "counts.assembly")

counts.gencode.GUSTO <- data.frame(counts(GUSTO.drgC.gencode.gdm[[3]]))
counts.gencode.GUSTO <- counts.gencode.GUSTO[row.names(counts.gencode.GUSTO) %in% GNAS, ]
counts.gencode.GUSTO <- counts.gencode.GUSTO %>%
  rownames_to_column(var = "row_name") %>%
  pivot_longer(
    cols = -row_name,
    names_to = "column_name",
    values_to = "value")
names(counts.gencode.GUSTO) <- c("transcript_id", "SampleID", "counts.gencode")

# Combine count data and add GDM status
counts.GUSTO.GNAS <- counts.assembly.GUSTO
counts.GUSTO.GNAS$counts.gencode <- counts.gencode.GUSTO$counts.gencode
counts.GUSTO.GNAS <- data.frame(counts.GUSTO.GNAS)

counts.GUSTO.GNAS$col <- "A"
counts.GUSTO.GNAS[counts.GUSTO.GNAS$transcript_id == "ENST00000476196.5", "col"] <- "B"

counts.GUSTO.GNAS$col <- factor(counts.GUSTO.GNAS$col, levels = c("A", "B"))
counts.GUSTO.GNAS <- counts.GUSTO.GNAS[order(counts.GUSTO.GNAS$col), ]
counts.GUSTO.GNAS <- left_join(counts.GUSTO.GNAS, cordat.GUSTO[, c("SampleID", "gdm")])

counts.GUSTO.GNAS$col <- as.character(counts.GUSTO.GNAS$col)
counts.GUSTO.GNAS$gdm <- as.character(counts.GUSTO.GNAS$gdm)
counts.GUSTO.GNAS[counts.GUSTO.GNAS$gdm == 0 & counts.GUSTO.GNAS$col == "B", "col"] <- "C"

counts.GUSTO.GNAS$col <- factor(counts.GUSTO.GNAS$col, levels = c("A", "C", "B"))
counts.GUSTO.GNAS <- counts.GUSTO.GNAS[order(counts.GUSTO.GNAS$col), ]

## Figure 3F: Read count correlation for GNAS
fig5f <- ggplot(counts.GUSTO.GNAS, aes(x = log2(counts.gencode + 1),
                                       y = log2(counts.assembly + 1),
                                       color = col)) +
  geom_point() + 
  theme_classic() +
  scale_color_manual(values = c("blue", "red", "grey50"),
                     limits = c("B", "C"),
                     name = "GDM",
                     labels = c("1", "0")) +
  xlab("GENCODEv45 counts (log2)") +
  ylab("lr-assembly counts (log2)") +
  ggtitle("Reads mapped to GNAS FSM\n(ENST00000476196.5 highlighted)") +
  theme(panel.grid.major.x = element_line(color = "grey95"),
        panel.grid.major.y = element_line(color = "grey95"),
        legend.position = "right") +
  geom_abline(slope = 1,
              intercept = 0,
              color = "black",
              linetype = "dashed") +
  xlim(0, 12) + ylim(0, 12)

fig5f
ggsave("Figures/figs5f.pdf", fig5f, width = 3.3, height = 3, units = "in", dpi = 300)

################################################################################
# Read Mapping Transition Analysis
################################################################################

## Parallel processing function for read mapping transitions
process_sample <- function(sample_id) {
  file_A <- file.path(data_dir, paste0(sample_id, "_gencode.tsv.gz"))   # gencode = Assembly A
  file_B <- file.path(data_dir, paste0(sample_id, "_assembly.tsv.gz"))  # assembly = Assembly B
  
  if (!file.exists(file_A) || !file.exists(file_B)) {
    message("Skipping missing files for: ", sample_id)
    return(NULL)
  }
  
  message("Processing sample: ", sample_id)
  
  message("Importing reads_A for: ", sample_id)
  reads_A <- as.data.table(read_tsv(file_A, col_names = c("read", "transcript_A"), show_col_types = FALSE))
  
  message("Importing reads_B for: ", sample_id)
  reads_B <- as.data.table(read_tsv(file_B, col_names = c("read", "transcript_B"), show_col_types = FALSE))
  
  message("Processing transitions for: ", sample_id)
  
  combined <- data.table(transcript_A = reads_A$transcript_A,
                         transcript_B = reads_B$transcript_B)
  
  rm(reads_A, reads_B)
  
  combined[, bucket_A := fifelse(
    transcript_A == "*", "unmapped",
    fifelse(transcript_A %in% target_A, transcript_A, "Other_A")
  )]
  
  setkey(combined, transcript_B)
  final <- dt_T[combined]
  rm(combined)
  
  final[is.na(category_B), category_B := "unmapped"]
  
  reads_summary <- data.frame(table(as.data.frame(final[, .(bucket_A, category_B)])))
  
  return(reads_summary)
}

# Set random seed for reproducible sample selection
set.seed(1234)
sample_ids_0 <- sample(cordat.GUSTO[cordat.GUSTO$gdm == 0, "SampleID"], 5)
sample_ids_1 <- sample(cordat.GUSTO[cordat.GUSTO$gdm == 1, "SampleID"], 5)

data_dir <- "DIR_ALIGN"

################################################################################
# Figure 5D: CSH1 Read Mapping Transitions
################################################################################

## Define CSH1 transcripts for transition analysis
target_A <- c("ENST00000453363.7", "ENST00000316193.13", "ENST00000329882.8",
              "ENST00000558284.1", "ENST00000610991.1", "ENST00000558661.1")

# Set up transcript categories for assembly
transcript_categories <- classif[, c("transcript_id", "structural_category")]
transcript_categories$structural_category <- as.character(transcript_categories$structural_category)
names(transcript_categories) <- c("transcript_B", "category_B")

target_B <- classif[classif$associated_gene == "ENSG00000136488.15", "transcript_id"]
transcript_categories[!transcript_categories$transcript_B %in% target_B, "category_B"] <- "non-CSH1"
transcript_categories[transcript_categories$transcript_B == "ENST00000558284.1", "category_B"] <- "ENST00000558284.1"
dt_T <- as.data.table(transcript_categories)
setkey(dt_T, transcript_B)

## Process GDM = 0 samples
flow_counts_all <- map_dfr(sample_ids_0, process_sample)

### Aggregate transitions across samples
flow_counts_combined.CSH1 <- flow_counts_all %>%
  group_by(bucket_A, category_B) %>%
  summarise(Freq = sum(Freq), .groups = "drop")

### Revalue factor levels for visualization
flow_counts_combined.CSH1$bucket_A <- revalue(flow_counts_combined.CSH1$bucket_A, 
                                              c("Other_A" = paste0("non-CSH1\n(N = ",
                                                                   length(tx2g.gencode$transcript_id) - length(target_A), ")")))

flow_counts_combined.CSH1$category_B <- revalue(flow_counts_combined.CSH1$category_B, 
                                                c("non-CSH1" = paste0("non-CSH1\n(N = ",
                                                                      length(transcript_categories[transcript_categories$category_B == "non-CSH1", "transcript_B"]), ")"),
                                                  "FSM" = paste0("FSM\n(N = ",
                                                                 length(transcript_categories[transcript_categories$category_B == "FSM", "transcript_B"]), ")"),
                                                  "ISM" = paste0("ISM\n(N = ",
                                                                 length(transcript_categories[transcript_categories$category_B == "ISM", "transcript_B"]), ")"),
                                                  "NIC" = paste0("NIC\n(N = ",
                                                                 length(transcript_categories[transcript_categories$category_B == "NIC", "transcript_B"]), ")"),
                                                  "NNC" = paste0("NNC\n(N = ",
                                                                 length(transcript_categories[transcript_categories$category_B == "NNC", "transcript_B"]), ")")))

flow_counts_combined.CSH1.GDM0 <- data.frame(flow_counts_combined.CSH1)

## Process GDM = 1 samples
flow_counts_all <- map_dfr(sample_ids_1, process_sample)

### Aggregate transitions across samples
flow_counts_combined.CSH1 <- flow_counts_all %>%
  group_by(bucket_A, category_B) %>%
  summarise(Freq = sum(Freq), .groups = "drop")

### Revalue factor levels
flow_counts_combined.CSH1$bucket_A <- revalue(flow_counts_combined.CSH1$bucket_A, 
                                              c("Other_A" = paste0("non-CSH1\n(N = ",
                                                                   length(tx2g.gencode$transcript_id) - length(target_A), ")")))

flow_counts_combined.CSH1$category_B <- revalue(flow_counts_combined.CSH1$category_B, 
                                                c("non-CSH1" = paste0("non-CSH1\n(N = ",
                                                                      length(transcript_categories[transcript_categories$category_B == "non-CSH1", "transcript_B"]), ")"),
                                                  "FSM" = paste0("FSM\n(N = ",
                                                                 length(transcript_categories[transcript_categories$category_B == "FSM", "transcript_B"]), ")"),
                                                  "ISM" = paste0("ISM\n(N = ",
                                                                 length(transcript_categories[transcript_categories$category_B == "ISM", "transcript_B"]), ")"),
                                                  "NIC" = paste0("NIC\n(N = ",
                                                                 length(transcript_categories[transcript_categories$category_B == "NIC", "transcript_B"]), ")"),
                                                  "NNC" = paste0("NNC\n(N = ",
                                                                 length(transcript_categories[transcript_categories$category_B == "NNC", "transcript_B"]), ")")))

flow_counts_combined.CSH1.GDM1 <- data.frame(flow_counts_combined.CSH1)

# Combine both GDM groups
flow_counts_combined.CSH1.GDM <- flow_counts_combined.CSH1.GDM0
flow_counts_combined.CSH1.GDM$Freq <- flow_counts_combined.CSH1.GDM$Freq + flow_counts_combined.CSH1.GDM1$Freq

# Prepare data for visualization (remove unmapped and low-frequency transitions)
flow_counts_combined.CSH1.GDM.plot <- flow_counts_combined.CSH1.GDM
flow_counts_combined.CSH1.GDM.plot <- flow_counts_combined.CSH1.GDM.plot[-c(37:48), ]
flow_counts_combined.CSH1.GDM.plot <- flow_counts_combined.CSH1.GDM.plot[-grep("unmapped", flow_counts_combined.CSH1.GDM.plot$category_B), ]

## Figure 5D: CSH1 read mapping transitions
fig5d <- ggplot(flow_counts_combined.CSH1.GDM.plot, 
                aes(axis1 = bucket_A, axis2 = category_B, y = sqrt(Freq))) +
  geom_alluvium(aes(fill = category_B), width = 1/4) +
  scale_fill_manual(values = c("red", myPalette[c(1, 3, 4)], "grey30")) +
  geom_stratum(width = 1/4, fill = "gray95", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), 
            size = 3, lineheight = 0.7) +
  scale_x_discrete(limits = c("gencode", "lr-assembly"), expand = c(.1, .1)) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "CSH1 Read Mapping Transitions", 
       y = expression(sqrt("Read Count")), x = NULL)

fig5d
ggsave("Figures/fig5d.pdf", fig5d, width = 8, height = 6, units = "in", dpi = 300, bg = 'white')

################################################################################
# Figure 5G: GNAS Read Mapping Transitions
################################################################################

## Define GNAS transcripts for transition analysis
target_A <- tx2g.gencode[tx2g.gencode$gene_id == "ENSG00000087460.29", "transcript_id"]

# Set up transcript categories for GNAS
transcript_categories <- classif[, c("transcript_id", "structural_category")]
transcript_categories$structural_category <- as.character(transcript_categories$structural_category)
names(transcript_categories) <- c("transcript_B", "category_B")

target_B <- classif[classif$associated_gene == "ENSG00000087460.29", "transcript_id"]
transcript_categories[!transcript_categories$transcript_B %in% target_B, "category_B"] <- "non-GNAS"
transcript_categories[transcript_categories$transcript_B == "ENST00000476196.5", "category_B"] <- "ENST00000476196.5"
dt_T <- as.data.table(transcript_categories)
setkey(dt_T, transcript_B)

## Process all samples for GNAS analysis
flow_counts_all <- map_dfr(c(sample_ids_0, sample_ids_1), process_sample)

## Aggregate transitions across samples
flow_counts_combined.GNAS <- flow_counts_all %>%
  group_by(bucket_A, category_B) %>%
  summarise(Freq = sum(Freq), .groups = "drop")

## Revalue factor levels for GNAS
flow_counts_combined.GNAS <- data.frame(flow_counts_combined.GNAS)
flow_counts_combined.GNAS$bucket_A <- as.character(flow_counts_combined.GNAS$bucket_A)

new_bucket_A <- target_A[!target_A == "ENST00000476196.5"]

flow_counts_combined.GNAS[flow_counts_combined.GNAS$bucket_A %in% new_bucket_A, "bucket_A"] <- "Other GNAS\n(N = 76)"
flow_counts_combined.GNAS$bucket_A <- factor(flow_counts_combined.GNAS$bucket_A)
flow_counts_combined.GNAS$bucket_A <- revalue(flow_counts_combined.GNAS$bucket_A, 
                                              c("Other_A" = paste0("non-GNAS\n(N = ",
                                                                   length(tx2g.gencode$transcript_id) - length(target_A), ")")))

flow_counts_combined.GNAS$category_B <- revalue(flow_counts_combined.GNAS$category_B, 
                                                c("non-GNAS" = paste0("non-GNAS\n(N = ",
                                                                      length(transcript_categories[transcript_categories$category_B == "non-GNAS", "transcript_B"]), ")"),
                                                  "FSM" = paste0("FSM\n(N = ",
                                                                 length(transcript_categories[transcript_categories$category_B == "FSM", "transcript_B"]), ")"),
                                                  "ISM" = paste0("ISM\n(N = ",
                                                                 length(transcript_categories[transcript_categories$category_B == "ISM", "transcript_B"]), ")"),
                                                  "NIC" = paste0("NIC\n(N = ",
                                                                 length(transcript_categories[transcript_categories$category_B == "NIC", "transcript_B"]), ")"),
                                                  "NNC" = paste0("NNC\n(N = ",
                                                                 length(transcript_categories[transcript_categories$category_B == "NNC", "transcript_B"]), ")")))

# Prepare data for visualization
flow_counts_combined.GNAS.plot <- flow_counts_combined.GNAS
flow_counts_combined.GNAS.plot <- flow_counts_combined.GNAS.plot[-c(309:316), ]
flow_counts_combined.GNAS.plot <- flow_counts_combined.GNAS.plot[-grep("unmapped", flow_counts_combined.GNAS.plot$category_B), ]

# Adjust frequencies for better visualization
flow_counts_combined.GNAS.plot[flow_counts_combined.GNAS.plot$bucket_A == "ENST00000476196.5", "Freq"] <- 
  flow_counts_combined.GNAS.plot[flow_counts_combined.GNAS.plot$bucket_A == "ENST00000476196.5", "Freq"] * 2
flow_counts_combined.GNAS.plot[grep("non-GNAS", flow_counts_combined.GNAS.plot$category_B), "Freq"] <- 
  flow_counts_combined.GNAS.plot[grep("non-GNAS", flow_counts_combined.GNAS.plot$category_B), "Freq"] * 8

flow_counts_combined.GNAS.plot0 <- flow_counts_combined.GNAS.plot

flow_counts_combined.GNAS.plot <- flow_counts_combined.GNAS.plot %>%
  dplyr::mutate(category_B = case_when(
    grepl("non-GNAS", category_B) ~ gsub("\n", "\n\n", as.character(category_B), fixed = TRUE),
    TRUE ~ as.character(category_B)
  ))

## Figure 5G: GNAS read mapping transitions
fig5g <- ggplot(flow_counts_combined.GNAS.plot, 
                aes(axis1 = bucket_A, axis2 = category_B, y = Freq)) +
  geom_alluvium(aes(fill = category_B), width = 1/4) +
  scale_fill_manual(values = c("red", myPalette[1], "grey30")) +
  geom_stratum(width = 1/4, fill = "gray95", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), 
            size = 3, lineheight = 0.7) +
  scale_x_discrete(limits = c("gencode", "lr-assembly"), expand = c(.1, .1)) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "GNAS Read Mapping Transitions", 
       y = expression(sqrt("Read Count")), x = NULL)

fig5g
ggsave("Figures/fig5g.pdf", fig5g, width = 8, height = 6, units = "in", dpi = 300, bg = 'white')


total_reads <- sum(flow_counts_combined.GNAS.plot$Freq)
non_gnas_reads <- sum(flow_counts_combined.GNAS.plot$Freq[grepl("non-GNAS", flow_counts_combined.GNAS.plot$category_B)])
percent_non_gnas <- (non_gnas_reads / total_reads) * 100

cat("Percent of reads mapping to non-GNAS:", round(percent_non_gnas, 2), "%\n")
cat("Non-GNAS reads:", non_gnas_reads, "\n")
cat("Total reads:", total_reads, "\n")

################################################################################
# End of Analysis
################################################################################