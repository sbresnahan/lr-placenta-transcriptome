################################################################################
# Differential Gene and Transcript Expression Analysis: GDM Association Study
# Author: Sean T. Bresnahan
# Description: Comprehensive differential expression analysis comparing gestational
#              diabetes mellitus (GDM) cases vs controls using both gene-level
#              (DEG) and transcript-level (DTE) approaches. Analyzes expression
#              patterns using filtered placental transcriptome assembly versus
#              GENCODE v45 reference across two independent cohorts.
#
# Input Files:
#   - Filtered assembly classification (ESPRESSO_corrected_SC_filtered_classification.txt)
#   - Filtered assembly GTF (ESPRESSO_corrected_SC_filtered.gtf)
#   - GENCODE v45 annotation GTF (gencode.v45.annotation.gtf)
#   - Salmon quantification results (assembly and GENCODE for both cohorts)
#   - GUSTO cohort metadata (20220823-Full_200_RNAseq_covars_v3.csv, v2.csv)
#   - Gen3G cohort metadata (multiple dbGaP files)
#   - OGTT results (rna_ogtt_results.xlsx)
#
# Output Data:
#   - GUSTO_ESPRESSO_assembly_filtered_SC.RData: GUSTO differential expression results
#   - Gen3G_ESPRESSO_assembly_filtered_SC.RData: Gen3G differential expression results
#
# Output Figures:
#   - Main Manuscript: 3A, 3B, 3C, 3D, 3E, 3F, 3G
#     * 3A: Venn diagrams comparing DTE/DGE between cohorts and annotations
#     * 3B: Log fold change comparison for CSH1 full-splice matches
#     * 3C: Read count correlation for CSH1 transcripts by GDM status
#     * 3D: Read mapping transitions for CSH1 between annotations
#     * 3E: Log fold change comparison for GNAS full-splice matches
#     * 3F: Read count correlation for GNAS transcripts by GDM status
#     * 3G: Read mapping transitions for GNAS between annotations
#
# Methods:
#   - RUVSeq batch correction with empirical control genes
#   - DESeq2 differential expression analysis
#   - Fishpond QU (quasi-UMI) correction for technical variance
#   - Mixed-effects modeling for covariates (sex, gestational age, batch)
#   - Isotwas gene-level aggregation for transcript-level results
#   - Alluvial flow analysis for read mapping transitions
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
    
    ttt <- merge(res, gene_res[, c('gene_id', 'Screen.P.Adj')], by = 'gene_id')
    
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
classif <- read.table("ESPRESSO_corrected_SC_filtered_classification.txt", 
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

load("tx2g.RData")

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

## Import and process GUSTO metadata
metadata <- read.csv("20220823-Full_200_RNAseq_covars_v3.csv")
metadata$SampleID <- paste0("J", metadata$ID)
metadata <- metadata[, c(51, 1:50)]
row.names(metadata) <- metadata$SampleID

# Add RNA quality metrics
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
names(cordat)[which(names(cordat) == "c_weight_birth_kg")] <- "birth_weight"

# Import OGTT results from supplementary analysis
gdm_cont <- readxl::read_excel('rna_ogtt_results.xlsx')
gdm_cont$SampleID <- paste0("J", gdm_cont$ID)
gdm_cont$ogtt_2hour_pw26 <- as.numeric(scale(gdm_cont$ogtt_2hour_pw26))

################################################################################
# GUSTO Assembly Quantification and Differential Expression
################################################################################

## Prepare assembly gene expression data
dirs.salmon.assembly <- list.dirs("COUNTS_assembly_SC", recursive = FALSE)
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
save(list=c("metadata.GUSTO","cordat.GUSTO"),file="metadata_GUSTO.RData")

### Gene-level differential expression analysis
drgC.assembly.gdm.gene <- deseq2_ruv(txi.assembly.filter, cordat, tx = FALSE, tx2g.assembly,
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
drgC.assembly.gdm <- deseq2_ruv(txi.assembly.filter, cordat, tx = TRUE, tx2g.assembly,
                                "gdm", "categorical", c("0", "1"), "SampleID",
                                "Batch", weights = TRUE, W = 2)

################################################################################
# GUSTO GENCODE Analysis
################################################################################

## Gene-level expression analysis with GENCODE
dirs.salmon.gencode <- list.dirs("COUNTS_gencode", recursive = FALSE)
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
drgC.gencode.gdm.gene <- deseq2_ruv(txi.gencode.filter, cordat, tx = FALSE, tx2g.gencode,
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
drgC.gencode.gdm <- deseq2_ruv(txi.gencode.filter, cordat, tx = TRUE, tx2g.gencode,
                               "gdm", "categorical", c("0", "1"), "SampleID",
                               "Batch", weights = TRUE, W = 2)

# Save GUSTO results for cross-cohort comparison
GUSTO.drgC.assembly.gdm <- drgC.assembly.gdm
GUSTO.drgC.assembly.gdm.gene <- drgC.assembly.gdm.gene
GUSTO.drgC.gencode.gdm <- drgC.gencode.gdm
GUSTO.drgC.gencode.gdm.gene <- drgC.gencode.gdm.gene
save(list = c("GUSTO.drgC.assembly.gdm",
              "GUSTO.drgC.gencode.gdm",
              "GUSTO.drgC.assembly.gdm.gene",
              "GUSTO.drgC.gencode.gdm.gene"),
     file = "GUSTO_ESPRESSO_assembly_filtered_SC.RData")

################################################################################
# Gen3G Cohort Analysis
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
# Gen3G Assembly Analysis
################################################################################

## Prepare assembly gene expression data
dirs.salmon.assembly <- list.dirs("COUNTS_assembly_SC", 
                                  recursive = FALSE)
files.list.assembly <- list()
for (i in 1:length(dirs.salmon.assembly)) {
  files.list.assembly[i] <- paste(unlist(dirs.salmon.assembly[i]), "quant.sf", sep = "/")
}
dirnames.assembly <- sapply(strsplit(unlist(dirs.salmon.assembly), "/"), "[",
                            length(unlist(strsplit(unlist(dirs.salmon.assembly), "/")[1])))
names(files.list.assembly) <- dirnames.assembly

# Filter to samples present in metadata
keep <- which(dirnames.assembly %in% cordat$Run)
dirs.salmon.assembly <- dirs.salmon.assembly[keep]
dirnames.assembly <- dirnames.assembly[keep]
files.list.assembly <- unlist(files.list.assembly)[keep]

# Import gene-level quantification
txi.assembly.g <- tximport(unlist(files.list.assembly), type = "salmon", tx2gene = tx2g.assembly,
                           geneIdCol = "gene_id", txIdCol = "transcript_id", dropInfReps = TRUE, txOut = FALSE)

### Quality control and batch filtering
dds.assembly <- DESeqDataSetFromTximport(txi.assembly.g, cordat, ~1)
vst.assembly <- vst(dds.assembly, blind = TRUE)
pca.assembly.check <- plotPCA(vst.assembly, intgroup = "sex")
pca.assembly.check

# Check batch effects
boxplot(gdm ~ sequencing_batch, cordat)
data.frame(table(cordat[, c("gdm", "sequencing_batch")]))

# Filter out problematic sequencing batches
filter.assembly <- cordat[!cordat$sequencing_batch %in% c("SK-3KKV", "SK-3KKW"), "Run"]
txi.assembly.filter <- filter_tximport(txi.assembly.g, filter.assembly)
cordat <- cordat[cordat$sequencing_batch %in% c("SK-3KKV", "SK-3KKW"), ]

# Save metadata for later use
metadata.Gen3G <- metadata
cordat.Gen3G <- cordat
save(list=c("metadata.Gen3G","cordat.Gen3G"),file="metadata_Gen3G.RData")

# Gene-level differential expression
drgC.assembly.gdm.gene <- deseq2_ruv(txi.assembly.filter, cordat, tx = FALSE, tx2g.assembly,
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
drgC.assembly.gdm <- deseq2_ruv(txi.assembly.filter, cordat, tx = TRUE, tx2g.assembly,
                                "gdm", "categorical", c("0", "1"), "Run",
                                "sequencing_batch", weights = TRUE, W = 2)

################################################################################
# Gen3G GENCODE Analysis
################################################################################

## Prepare GENCODE gene expression data
dirs.salmon.gencode <- list.dirs("COUNTS_gencode", 
                                 recursive = FALSE)
files.list.gencode <- list()
for (i in 1:length(dirs.salmon.gencode)) {
  files.list.gencode[i] <- paste(unlist(dirs.salmon.gencode[i]), "quant.sf", sep = "/")
}
dirnames.gencode <- sapply(strsplit(unlist(dirs.salmon.gencode), "/"), "[",
                           length(unlist(strsplit(unlist(dirs.salmon.gencode), "/")[1])))
names(files.list.gencode) <- dirnames.gencode

# Filter to samples present in metadata
keep <- which(dirnames.gencode %in% cordat$Run)
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
drgC.gencode.gdm.gene <- deseq2_ruv(txi.gencode.filter, cordat, tx = FALSE, tx2g.gencode,
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
drgC.gencode.gdm <- deseq2_ruv(txi.gencode.filter, cordat, tx = TRUE, tx2g.gencode,
                               "gdm", "categorical", c("0", "1"), "Run",
                               "sequencing_batch", weights = TRUE, W = 2)

# Save Gen3G results for cross-cohort comparison
Gen3G.drgC.assembly.gdm <- drgC.assembly.gdm
Gen3G.drgC.assembly.gdm.gene <- drgC.assembly.gdm.gene
Gen3G.drgC.gencode.gdm <- drgC.gencode.gdm
Gen3G.drgC.gencode.gdm.gene <- drgC.gencode.gdm.gene
save(list = c("Gen3G.drgC.assembly.gdm", "Gen3G.drgC.gencode.gdm",
              "Gen3G.drgC.assembly.gdm.gene", "Gen3G.drgC.gencode.gdm.gene"),
     file = "Gen3G_ESPRESSO_assembly_filtered_SC.RData")

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
Gen3G.vst.assembly <- assay(varianceStabilizingTransformation(Gen3G.drgC.assembly.gdm[[3]], blind = FALSE))
Gen3G.vst.gencode <- assay(varianceStabilizingTransformation(Gen3G.drgC.gencode.gdm[[3]], blind = FALSE))

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

results.assembly.GUSTO.gdm.g <- GUSTO.drgC.assembly.gdm[[5]]
results.assembly.Gen3G.gdm.g <- Gen3G.drgC.assembly.gdm[[5]]
results.gencode.GUSTO.gdm.g <- GUSTO.drgC.gencode.gdm[[5]]
results.gencode.Gen3G.gdm.g <- Gen3G.drgC.gencode.gdm[[5]]

results.assembly.GUSTO.gdm.gene <- GUSTO.drgC.assembly.gdm.gene[[4]]
results.assembly.GUSTO.gdm.gene$gene_id <- row.names(results.assembly.GUSTO.gdm.gene)
results.assembly.Gen3G.gdm.gene <- Gen3G.drgC.assembly.gdm.gene[[4]]
results.assembly.Gen3G.gdm.gene$gene_id <- row.names(results.assembly.Gen3G.gdm.gene)

results.gencode.GUSTO.gdm.gene <- GUSTO.drgC.gencode.gdm.gene[[4]]
results.gencode.GUSTO.gdm.gene$gene_id <- row.names(results.gencode.GUSTO.gdm.gene)
results.gencode.Gen3G.gdm.gene <- Gen3G.drgC.gencode.gdm.gene[[4]]
results.gencode.Gen3G.gdm.gene$gene_id <- row.names(results.gencode.Gen3G.gdm.gene)

################################################################################
# Significant Gene/Transcript Lists
################################################################################

## Total numbers of DTEs and DGEs at FDR < 0.1
GUSTO.gdm.DTE.assembly <- results.assembly.GUSTO.gdm[results.assembly.GUSTO.gdm$padj < 0.1, "transcript_id"]
GUSTO.gdm.DTE.gencode <- results.gencode.GUSTO.gdm[results.gencode.GUSTO.gdm$padj < 0.1, "transcript_id"]
Gen3G.gdm.DTE.assembly <- results.assembly.Gen3G.gdm[results.assembly.Gen3G.gdm$padj < 0.1, "transcript_id"]
Gen3G.gdm.DTE.gencode <- results.gencode.Gen3G.gdm[results.gencode.Gen3G.gdm$padj < 0.1, "transcript_id"]

GUSTO.gdm.DGE.assembly <- results.assembly.GUSTO.gdm.gene[results.assembly.GUSTO.gdm.gene$padj < 0.1, "gene_id"]
GUSTO.gdm.DGE.gencode <- results.gencode.GUSTO.gdm.gene[results.gencode.GUSTO.gdm.gene$padj < 0.1, "gene_id"]
Gen3G.gdm.DGE.assembly <- results.assembly.Gen3G.gdm.gene[results.assembly.Gen3G.gdm.gene$padj < 0.1, "gene_id"]
Gen3G.gdm.DGE.gencode <- results.gencode.Gen3G.gdm.gene[results.gencode.Gen3G.gdm.gene$padj < 0.1, "gene_id"]

################################################################################
# Figure 3A: Venn Diagrams Comparing DTE and DGE
################################################################################

# DTE comparison across cohorts and annotations
venn_list <- list(
  "GUSTO\nlr-assembly\n" = GUSTO.gdm.DTE.assembly,
  "GUSTO\ngencode\n" = GUSTO.gdm.DTE.gencode,
  "Gen3G\nlr-assembly\n" = Gen3G.gdm.DTE.assembly,
  "Gen3G\ngencode\n" = Gen3G.gdm.DTE.gencode
)

venn_list <- lapply(venn_list, as.character)

## Figure 3A left panel: DTE comparison
fig3a1 <- ggVennDiagram(venn_list, label_alpha = 0, label = "count") +
  theme_void() +
  scale_fill_gradient(low = "#FFFFFF", high = "grey") +
  ggtitle("DTE") + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(expand = expansion(mult = .2)) +
  scale_x_continuous(expand = expansion(mult = .2))

fig3a1
ggsave("fig3a1.png", fig3a1, width = 5, height = 5, units = "in", dpi = 300)

### DGE comparison across cohorts and annotations
venn_list <- list(
  "GUSTO\nlr-assembly\n" = GUSTO.gdm.DGE.assembly,
  "GUSTO\ngencode\n" = GUSTO.gdm.DGE.gencode,
  "Gen3G\nlr-assembly\n" = Gen3G.gdm.DGE.assembly,
  "Gen3G\ngencode\n" = Gen3G.gdm.DGE.gencode
)

venn_list <- lapply(venn_list, as.character)

## Figure 3A right panel: DGE comparison
fig3a2 <- ggVennDiagram(venn_list, label_alpha = 0, label = "count") +
  theme_void() +
  scale_fill_gradient(low = "#FFFFFF", high = "grey") +
  ggtitle("DGE") + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(expand = expansion(mult = .2)) +
  scale_x_continuous(expand = expansion(mult = .2))

fig3a2
ggsave("fig3a2.png", fig3a2, width = 5, height = 5, units = "in", dpi = 300)

################################################################################
# Figure 3B: CSH1 Log Fold Change Comparison
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
fig3b <- ggplot(dat.lfc.csh1, aes(x = lfc.gencode, y = lfc.assembly)) +
  theme_classic() +
  xlab("logFC (gencode)") +
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

fig3b
ggsave("fig3b.png", fig3b, width = 3, height = 3, units = "in", dpi = 300)

################################################################################
# Figure 3C: CSH1 Read Count Correlation
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
fig3c <- ggplot(counts.GUSTO.CSH1, aes(x = log2(counts.gencode + 1),
                                       y = log2(counts.assembly + 1),
                                       color = col)) +
  geom_point() + 
  theme_classic() +
  scale_color_manual(values = c("blue", "red", "grey50"),
                     limits = c("B", "C"),
                     name = "GDM",
                     labels = c("1", "0")) +
  xlab("gencode counts (log2)") +
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

fig3c
ggsave("fig3c.png", fig3c, width = 3.5, height = 3, units = "in", dpi = 300)

################################################################################
# Figure 3E: GNAS Log Fold Change Comparison
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
fig3e <- ggplot(dat.lfc.GNAS, aes(x = lfc.gencode, y = lfc.assembly)) +
  theme_classic() +
  xlab("logFC (gencode)") +
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

fig3e
ggsave("fig3e.png", fig3e, width = 3, height = 3, units = "in", dpi = 300)

################################################################################
# Figure 3F: GNAS Read Count Correlation
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
fig3f <- ggplot(counts.GUSTO.GNAS, aes(x = log2(counts.gencode + 1),
                                       y = log2(counts.assembly + 1),
                                       color = col)) +
  geom_point() + 
  theme_classic() +
  scale_color_manual(values = c("blue", "red", "grey50"),
                     limits = c("B", "C"),
                     name = "GDM",
                     labels = c("1", "0")) +
  xlab("gencode counts (log2)") +
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

fig3f
ggsave("fig3f.png", fig3f, width = 3.5, height = 3, units = "in", dpi = 300)

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
# Figure 3D: CSH1 Read Mapping Transitions
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

## Figure 3D: CSH1 read mapping transitions
fig3d <- ggplot(flow_counts_combined.CSH1.GDM.plot, 
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

fig3d
ggsave("fig3d.png", fig3d, width = 7.74, height = 6.17, units = "in", dpi = 300, bg = 'white')

################################################################################
# Figure 3G: GNAS Read Mapping Transitions
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

## Figure 3G: GNAS read mapping transitions
fig3g <- ggplot(flow_counts_combined.GNAS.plot, 
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

fig3g
ggsave("fig3g.png", fig3g, width = 7.74, height = 6.17, units = "in", dpi = 300, bg = 'white')

################################################################################
# End of Analysis
################################################################################