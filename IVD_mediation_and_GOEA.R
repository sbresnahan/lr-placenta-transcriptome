################################################################################
# Mediation Analysis and Gene Ontology Enrichment Analysis: GDM Association Study
# Author: Sean T. Bresnahan
# Description: Mediation analysis investigating how differentially expressed 
#              transcripts mediate the relationship between gestational diabetes 
#              mellitus (GDM) and birth weight. Includes individual transcript 
#              mediation analysis and principal component-based mediation analysis.
#              Also performs Gene Ontology enrichment analysis on significant 
#              transcript sets using enrichR.
#
# Input Files:
#   - tx2g.RData: Transcript-to-gene mapping tables
#   - GUSTO_ESPRESSO_assembly_filtered_SC.RData: GUSTO differential expression results
#   - Gen3G_ESPRESSO_assembly_filtered_SC.RData: Gen3G differential expression results
#   - ESPRESSO_corrected_SC_filtered_classification.txt: Assembly structural classifications
#   - Metadata files for both cohorts (same as IVC_DGE_DTE.R)
#   - ESPRESSO_corrected_SC_filtered.gtf: Assembly annotation
#   - gencode.v45.annotation.gtf: Reference annotation
#
# Output Data:
#   - Mediation analysis results for individual transcripts
#   - Principal component mediation analysis results
#   - Gene Ontology enrichment analysis results (CSV files)
#
# Output Figures:
#   - Figure 4B: PC mediation sensitivity analysis
#   - Figure 4C: Top transcript mediators (GUSTO lr-assembly)
#   - Figure 4D: Top transcript mediators (Gen3G lr-assembly)
#   - Figure 4E: CSH1 transcript isoforms mediating GDM effects
#   - Figure S8: PC mediation sensitivity analysis for Gen3G
#   - EnrichR lollipop plots for pathway enrichment
#
# Methods:
#   - Individual transcript mediation using 'mediation' package
#   - Principal component-based mediation using structural equation modeling (lavaan)
#   - Gene Ontology enrichment analysis using enrichR
#   - Visualization of transcript structures using ggtranscript
################################################################################

################################################################################
# Environment Setup
################################################################################

# Configure R library paths for HPC environment
.libPaths(c("/home/stbresnahan/R/ubuntu/4.3.1", .libPaths()))

# Load required libraries
library(dplyr)           # Data manipulation
library(rtracklayer)     # Genomic data import/export
library(ggplot2)         # Data visualization
library(tidyr)           # Data tidying
library(cowplot)         # Plot arrangement
library(plyr)            # Data manipulation (legacy functions)
library(scales)          # Scale functions for ggplot2
library(ggExtra)         # Additional ggplot2 functionality
library(edgeR)           # RNA-seq analysis toolkit
library(tximport)        # Transcript abundance import
library(DESeq2)          # Differential expression analysis
library(fishpond)        # Transcript-level inferential variance
library(tximeta)         # Metadata-aware transcript quantification import
library(ggdist)          # Distribution visualization
library(ggpubr)          # Publication-ready plots
library(tibble)          # Enhanced data frames
library(mediation)       # Mediation analysis
library(foreach)         # Parallel iteration
library(doParallel)      # Parallel processing backend
library(parallel)        # Parallel processing utilities
library(lavaan)          # Structural equation modeling
library(effectsize)      # Effect size calculations
library(data.table)      # High-performance data manipulation
library(enrichR)         # Gene enrichment analysis
library(org.Hs.eg.db)    # Human genome annotation
library(forcats)         # Factor manipulation
library(ggbreak)         # Axis break functionality
library(GO.db)           # Gene Ontology database
library(AnnotationDbi)   # Annotation database interface
library(biomaRt)         # Biomart interface
library(ggtranscript)    # Transcript structure visualization

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

# Load transcript-to-gene mapping tables
load("tx2g.RData")

# Load short-read cohort metadata
load("metadata_GUSTO.RData")
load("metadata_Gen3G.RData")

# Load differential expression results from both cohorts
load("GUSTO_ESPRESSO_assembly_filtered_SC.RData")
load("Gen3G_ESPRESSO_assembly_filtered_SC.RData")

# Import assembly isoform classification from SQANTI3
classif <- read.table("ESPRESSO_corrected_SC_filtered_classification.txt", 
                      header = TRUE)

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

names(classif)[1] <- "transcript_id"

# Create simplified classification for pathway analysis
classif.pcorr <- classif
classif.pcorr$structural_category <- revalue(classif.pcorr$structural_category, 
                                             c("Genic Genomic" = "Other",
                                               "Antisense" = "Other",
                                               "Fusion" = "Other",
                                               "Intergenic" = "Other",
                                               "Genic Intron" = "Other"))

################################################################################
# Variance Stabilized Transformation for Expression Data
################################################################################

# Generate VST transformed read counts for visualization and mediation analysis
GUSTO.vst.assembly <- assay(varianceStabilizingTransformation(GUSTO.drgC.assembly.gdm[[3]], blind = FALSE))
GUSTO.vst.gencode <- assay(varianceStabilizingTransformation(GUSTO.drgC.gencode.gdm[[3]], blind = FALSE))
Gen3G.vst.assembly <- assay(varianceStabilizingTransformation(Gen3G.drgC.assembly.gdm[[3]], blind = FALSE))
Gen3G.vst.gencode <- assay(varianceStabilizingTransformation(Gen3G.drgC.gencode.gdm[[3]], blind = FALSE))

################################################################################
# Results Tables Preparation
################################################################################

# Extract differential expression results for transcript-level analysis
results.assembly.GUSTO.gdm <- GUSTO.drgC.assembly.gdm[[4]]
results.assembly.GUSTO.gdm$transcript_id <- row.names(results.assembly.GUSTO.gdm)
results.assembly.Gen3G.gdm <- Gen3G.drgC.assembly.gdm[[4]]
results.assembly.Gen3G.gdm$transcript_id <- row.names(results.assembly.Gen3G.gdm)

results.gencode.GUSTO.gdm <- GUSTO.drgC.gencode.gdm[[4]]
results.gencode.GUSTO.gdm$transcript_id <- row.names(results.gencode.GUSTO.gdm)
results.gencode.Gen3G.gdm <- Gen3G.drgC.gencode.gdm[[4]]
results.gencode.Gen3G.gdm$transcript_id <- row.names(results.gencode.Gen3G.gdm)

# Extract gene-level aggregated results for transcript-level analysis
results.assembly.GUSTO.gdm.g <- GUSTO.drgC.assembly.gdm[[5]]
results.assembly.Gen3G.gdm.g <- Gen3G.drgC.assembly.gdm[[5]]
results.gencode.GUSTO.gdm.g <- GUSTO.drgC.gencode.gdm[[5]]
results.gencode.Gen3G.gdm.g <- Gen3G.drgC.gencode.gdm[[5]]

# Extract gene-level differential expression results
results.assembly.GUSTO.gdm.gene <- GUSTO.drgC.assembly.gdm.gene[[4]]
results.assembly.GUSTO.gdm.gene$gene_id <- row.names(results.assembly.GUSTO.gdm.gene)
results.assembly.Gen3G.gdm.gene <- Gen3G.drgC.assembly.gdm.gene[[4]]
results.assembly.Gen3G.gdm.gene$gene_id <- row.names(results.assembly.Gen3G.gdm.gene)

results.gencode.GUSTO.gdm.gene <- GUSTO.drgC.gencode.gdm.gene[[4]]
results.gencode.GUSTO.gdm.gene$gene_id <- row.names(results.gencode.GUSTO.gdm.gene)
results.gencode.Gen3G.gdm.gene <- Gen3G.drgC.gencode.gdm.gene[[4]]
results.gencode.Gen3G.gdm.gene$gene_id <- row.names(results.gencode.Gen3G.gdm.gene)

# Define significant gene/transcript lists at FDR < 0.1
GUSTO.gdm.DTE.assembly <- results.assembly.GUSTO.gdm[results.assembly.GUSTO.gdm$padj < 0.1, "transcript_id"]
GUSTO.gdm.DTE.gencode <- results.gencode.GUSTO.gdm[results.gencode.GUSTO.gdm$padj < 0.1, "transcript_id"]
Gen3G.gdm.DTE.assembly <- results.assembly.Gen3G.gdm[results.assembly.Gen3G.gdm$padj < 0.1, "transcript_id"]
Gen3G.gdm.DTE.gencode <- results.gencode.Gen3G.gdm[results.gencode.Gen3G.gdm$padj < 0.1, "transcript_id"]

GUSTO.gdm.DGEtx.assembly <- unique(data.frame(results.assembly.GUSTO.gdm.g[results.assembly.GUSTO.gdm.g$Screen.P.Adj < 0.1, "gene_id"])[, 1])
GUSTO.gdm.DGEtx.gencode <- unique(data.frame(results.gencode.GUSTO.gdm.g[results.gencode.GUSTO.gdm.g$Screen.P.Adj < 0.1, "gene_id"])[, 1])
Gen3G.gdm.DGEtx.assembly <- unique(data.frame(results.assembly.Gen3G.gdm.g[results.assembly.Gen3G.gdm.g$Screen.P.Adj < 0.1, "gene_id"])[, 1])
Gen3G.gdm.DGEtx.gencode <- unique(data.frame(results.gencode.Gen3G.gdm.g[results.gencode.Gen3G.gdm.g$Screen.P.Adj < 0.1, "gene_id"])[, 1])

GUSTO.gdm.DGE.assembly <- results.assembly.GUSTO.gdm.gene[results.assembly.GUSTO.gdm.gene$padj < 0.1, "gene_id"]
GUSTO.gdm.DGE.gencode <- results.gencode.GUSTO.gdm.gene[results.gencode.GUSTO.gdm.gene$padj < 0.1, "gene_id"]
Gen3G.gdm.DGE.assembly <- results.assembly.Gen3G.gdm.gene[results.assembly.Gen3G.gdm.gene$padj < 0.1, "gene_id"]
Gen3G.gdm.DGE.gencode <- results.gencode.Gen3G.gdm.gene[results.gencode.Gen3G.gdm.gene$padj < 0.1, "gene_id"]

################################################################################
# Individual Transcript Mediation Analysis
################################################################################

# Set up parallel backend for computationally intensive mediation analysis
num_cores <- detectCores() - 1  # Leave one core free for system processes
cl <- makeCluster(num_cores)
registerDoParallel(cl)

## GUSTO assembly: Individual transcript mediation analysis
GUSTO.expr.assembly <- data.frame(GUSTO.vst.assembly)
GUSTO.expr.assembly$transcript_id <- row.names(GUSTO.expr.assembly)

# Filter to significant differentially expressed transcripts
GUSTO.expr.assembly <- GUSTO.expr.assembly[GUSTO.expr.assembly$transcript_id %in% 
                                             results.assembly.GUSTO.gdm[results.assembly.GUSTO.gdm$padj < 0.1, "transcript_id"], ]

# Reshape data for mediation analysis
GUSTO.expr.assembly <- GUSTO.expr.assembly %>%
  pivot_longer(
    cols = -transcript_id,
    names_to = "SampleID",
    values_to = "TXE")

# Add covariate data
GUSTO.expr.assembly <- left_join(GUSTO.expr.assembly, 
                                 cordat.GUSTO[, c("SampleID", "gdm", "birth_weight", "sex", "GA")])
txids <- unique(GUSTO.expr.assembly$transcript_id)

# Parallel mediation analysis for each transcript
GUSTO.mediation.assembly <- foreach(tx = 1:length(txids), 
                                    .packages = c("mediation"),
                                    .combine = rbind) %dopar% {
                                      
                                      data <- GUSTO.expr.assembly[GUSTO.expr.assembly$transcript_id == txids[tx], ]
                                      mediator_model <- lm(TXE ~ gdm + sex + GA, data = data)
                                      outcome_model <- lm(birth_weight ~ gdm + TXE + sex + GA, data = data)
                                      
                                      tryCatch({
                                        mediate_result <- mediate(mediator_model, outcome_model, 
                                                                  treat = "gdm", mediator = "TXE", 
                                                                  boot = TRUE, sims = 10)
                                        data.frame(transcript_id = txids[tx],
                                                   estimate = mediate_result$d0,
                                                   lower_95CI = mediate_result$d0.ci[1],
                                                   upper_95CI = mediate_result$d0.ci[2],
                                                   pvalue = mediate_result$d0.p)
                                      }, error = function(e) {
                                        warning(paste("Error in transcript", txids[tx], ":", e$message))
                                        data.frame(transcript_id = txids[tx],
                                                   estimate = NA,
                                                   lower_95CI = NA,
                                                   upper_95CI = NA,
                                                   pvalue = NA)
                                      })
                                    }

row.names(GUSTO.mediation.assembly) <- NULL
GUSTO.mediation.assembly$estimate_abs <- abs(GUSTO.mediation.assembly$estimate)
GUSTO.mediation.assembly <- left_join(GUSTO.mediation.assembly, tx2g.assembly)

## GUSTO gencode: Individual transcript mediation analysis
GUSTO.expr.gencode <- data.frame(GUSTO.vst.gencode)
GUSTO.expr.gencode$transcript_id <- row.names(GUSTO.expr.gencode)

# Filter to significant differentially expressed transcripts
GUSTO.expr.gencode <- GUSTO.expr.gencode[GUSTO.expr.gencode$transcript_id %in% 
                                           results.gencode.GUSTO.gdm[results.gencode.GUSTO.gdm$padj < 0.1, "transcript_id"], ]

# Reshape data for mediation analysis
GUSTO.expr.gencode <- GUSTO.expr.gencode %>%
  pivot_longer(
    cols = -transcript_id,
    names_to = "SampleID",
    values_to = "TXE")

# Add covariate data
GUSTO.expr.gencode <- left_join(GUSTO.expr.gencode, 
                                cordat.GUSTO[, c("SampleID", "gdm", "birth_weight", "sex", "GA")])
txids <- unique(GUSTO.expr.gencode$transcript_id)

# Parallel mediation analysis for each transcript
GUSTO.mediation.gencode <- foreach(tx = 1:length(txids), 
                                   .packages = c("mediation"),
                                   .combine = rbind) %dopar% {
                                     
                                     data <- GUSTO.expr.gencode[GUSTO.expr.gencode$transcript_id == txids[tx], ]
                                     mediator_model <- lm(TXE ~ gdm + sex + GA, data = data)
                                     outcome_model <- lm(birth_weight ~ gdm + TXE + sex + GA, data = data)
                                     
                                     tryCatch({
                                       mediate_result <- mediate(mediator_model, outcome_model, 
                                                                 treat = "gdm", mediator = "TXE", 
                                                                 boot = TRUE, sims = 10)
                                       data.frame(transcript_id = txids[tx],
                                                  estimate = mediate_result$d0,
                                                  lower_95CI = mediate_result$d0.ci[1],
                                                  upper_95CI = mediate_result$d0.ci[2],
                                                  pvalue = mediate_result$d0.p)
                                     }, error = function(e) {
                                       warning(paste("Error in transcript", txids[tx], ":", e$message))
                                       data.frame(transcript_id = txids[tx],
                                                  estimate = NA,
                                                  lower_95CI = NA,
                                                  upper_95CI = NA,
                                                  pvalue = NA)
                                     })
                                   }

row.names(GUSTO.mediation.gencode) <- NULL
GUSTO.mediation.gencode$estimate_abs <- abs(GUSTO.mediation.gencode$estimate)
GUSTO.mediation.gencode <- left_join(GUSTO.mediation.gencode, tx2g.gencode)

# Combine GUSTO mediation results
GUSTO.mediation.assembly$annotation <- "lr-assembly"
GUSTO.mediation.gencode$annotation <- "gencode"

GUSTO.mediation.gdm <- rbind(GUSTO.mediation.gencode, GUSTO.mediation.assembly)
GUSTO.mediation.gdm[GUSTO.mediation.gdm$pvalue == 0, "pvalue"] <- 2e-16
GUSTO.mediation.gdm$cohort <- "GUSTO"

## Gen3G assembly: Individual transcript mediation analysis
Gen3G.expr.assembly <- data.frame(Gen3G.vst.assembly)
Gen3G.expr.assembly$transcript_id <- row.names(Gen3G.expr.assembly)

# Filter to significant differentially expressed transcripts
Gen3G.expr.assembly <- Gen3G.expr.assembly[Gen3G.expr.assembly$transcript_id %in% 
                                             results.assembly.Gen3G.gdm[results.assembly.Gen3G.gdm$padj < 0.1, "transcript_id"], ]

# Reshape data for mediation analysis
Gen3G.expr.assembly <- Gen3G.expr.assembly %>%
  pivot_longer(
    cols = -transcript_id,
    names_to = "Run",
    values_to = "TXE")

# Add covariate data (simpler model for Gen3G due to sample size)
Gen3G.expr.assembly <- left_join(Gen3G.expr.assembly, 
                                 cordat.Gen3G[, c("Run", "gdm", "birth_weight", "GA")])
txids <- unique(Gen3G.expr.assembly$transcript_id)

# Parallel mediation analysis for each transcript
Gen3G.mediation.assembly <- foreach(tx = 1:length(txids), 
                                    .packages = c("mediation"),
                                    .combine = rbind) %dopar% {
                                      
                                      data <- Gen3G.expr.assembly[Gen3G.expr.assembly$transcript_id == txids[tx], ]
                                      mediator_model <- lm(TXE ~ gdm, data = data)
                                      outcome_model <- lm(birth_weight ~ gdm + TXE, data = data)
                                      
                                      tryCatch({
                                        mediate_result <- mediate(mediator_model, outcome_model, treat = "gdm", mediator = "TXE", 
                                                                  boot = TRUE, sims = 100)
                                        data.frame(transcript_id = txids[tx],
                                                   estimate = mediate_result$d0,
                                                   lower_95CI = mediate_result$d0.ci[1],
                                                   upper_95CI = mediate_result$d0.ci[2],
                                                   pvalue = mediate_result$d0.p)
                                      }, error = function(e) {
                                        warning(paste("Error in transcript", txids[tx], ":", e$message))
                                        data.frame(transcript_id = txids[tx],
                                                   estimate = NA,
                                                   lower_95CI = NA,
                                                   upper_95CI = NA,
                                                   pvalue = NA)
                                      })
                                    }

row.names(Gen3G.mediation.assembly) <- NULL
Gen3G.mediation.assembly$estimate_abs <- abs(Gen3G.mediation.assembly$estimate)
Gen3G.mediation.assembly <- left_join(Gen3G.mediation.assembly, tx2g.assembly, by = "transcript_id")

## Gen3G gencode: Individual transcript mediation analysis
Gen3G.expr.gencode <- data.frame(Gen3G.vst.gencode)
Gen3G.expr.gencode$transcript_id <- row.names(Gen3G.expr.gencode)

# Filter to significant differentially expressed transcripts
Gen3G.expr.gencode <- Gen3G.expr.gencode[Gen3G.expr.gencode$transcript_id %in% 
                                           results.gencode.Gen3G.gdm[results.gencode.Gen3G.gdm$padj < 0.1, "transcript_id"], ]

# Reshape data for mediation analysis
Gen3G.expr.gencode <- Gen3G.expr.gencode %>%
  pivot_longer(
    cols = -transcript_id,
    names_to = "Run",
    values_to = "TXE")

# Add covariate data
Gen3G.expr.gencode <- left_join(Gen3G.expr.gencode, 
                                cordat.Gen3G[, c("Run", "gdm", "birth_weight", "GA")])
txids <- unique(Gen3G.expr.gencode$transcript_id)

# Parallel mediation analysis for each transcript
Gen3G.mediation.gencode <- foreach(tx = 1:length(txids), 
                                   .packages = c("mediation"),
                                   .combine = rbind) %dopar% {
                                     
                                     data <- Gen3G.expr.gencode[Gen3G.expr.gencode$transcript_id == txids[tx], ]
                                     mediator_model <- lm(TXE ~ gdm, data = data)
                                     outcome_model <- lm(birth_weight ~ gdm + TXE, data = data)
                                     
                                     tryCatch({
                                       mediate_result <- mediate(mediator_model, outcome_model, 
                                                                 treat = "gdm", mediator = "TXE", 
                                                                 boot = TRUE, sims = 100)
                                       data.frame(transcript_id = txids[tx],
                                                  estimate = mediate_result$d0,
                                                  lower_95CI = mediate_result$d0.ci[1],
                                                  upper_95CI = mediate_result$d0.ci[2],
                                                  pvalue = mediate_result$d0.p)
                                     }, error = function(e) {
                                       warning(paste("Error in transcript", txids[tx], ":", e$message))
                                       data.frame(transcript_id = txids[tx],
                                                  estimate = NA,
                                                  lower_95CI = NA,
                                                  upper_95CI = NA,
                                                  pvalue = NA)
                                     })
                                   }

row.names(Gen3G.mediation.gencode) <- NULL
Gen3G.mediation.gencode$estimate_abs <- abs(Gen3G.mediation.gencode$estimate)
Gen3G.mediation.gencode <- left_join(Gen3G.mediation.gencode, tx2g.gencode, by = "transcript_id")

# Combine Gen3G mediation results
Gen3G.mediation.assembly$annotation <- "lr-assembly"
Gen3G.mediation.gencode$annotation <- "gencode"

Gen3G.mediation.gdm <- rbind(Gen3G.mediation.gencode, Gen3G.mediation.assembly)
Gen3G.mediation.gdm[Gen3G.mediation.gdm$pvalue == 0, "pvalue"] <- 2e-16
Gen3G.mediation.gdm$cohort <- "Gen3G"

################################################################################
# Principal Component-Based Mediation Analysis Functions
################################################################################

# Function to perform PC sensitivity analysis across different numbers of components
PCmed.sensitivity <- function(expmat, cordat, PCmax) {
  
  require(lavaan)
  require(dplyr)
  
  # Helper function to generate mediation model syntax
  generate_mediation_model_body <- function(n_pcs) {
    mediator_block <- paste0(
      "  PC", 1:n_pcs, " ~ a", 1:n_pcs, "*gdm + sex_m", 1:n_pcs, "*sex + GA_m", 1:n_pcs, "*GA + ethnicity_m", 1:n_pcs, "*mother_ethnicity",
      collapse = "\n"
    )
    
    outcome_block <- paste0(
      "  birth_weight ~ ", paste0("b", 1:n_pcs, "*PC", 1:n_pcs, collapse = " + ")
    )
    
    indirect_block <- paste0(
      "  ind", 1:n_pcs, " := a", 1:n_pcs, "*b", 1:n_pcs, collapse = "\n"
    )
    
    total_indirect_block <- paste0(
      "  total_indirect := ", paste0("ind", 1:n_pcs, collapse = " + ")
    )
    
    total_effect_block <- "  total := c + total_indirect"
    prop_mediated_block <- "  prop_mediated := total_indirect/total"
    
    prop_each_block <- paste0(
      "  prop_med", 1:n_pcs, " := ind", 1:n_pcs, "/total", collapse = "\n"
    )
    
    paste(
      "# Mediator paths with sex, GA, and ethnicity adjustment\n", mediator_block,
      "\n\n# Outcome paths\n", outcome_block,
      "\n\n# Compute indirect effects\n", indirect_block,
      "\n\n# Total indirect effect\n", total_indirect_block,
      "\n\n# Total effect\n", total_effect_block,
      "\n\n# Proportion mediated (overall)\n", prop_mediated_block,
      "\n\n# Proportion for each mediator\n", prop_each_block,
      "\n", sep = ""
    )
  }
  
  # Perform PCA and extract variance explained
  pca_res <- prcomp(t(as.matrix(expmat)))
  pca_df <- as.data.frame(pca_res$x[, 1:PCmax])  # top PCmax PCs
  pca_var_explained <- summary(pca_res)$importance["Proportion of Variance", 1:PCmax]
  pca_cumvar <- cumsum(pca_var_explained)
  
  results <- list()
  
  # Test different numbers of principal components
  for (i in 1:PCmax) {
    # Prepare data with top i PCs
    df <- data.frame(pca_df[, 1:i])
    colnames(df) <- paste0("PC", 1:i)
    df$SampleID <- rownames(pca_df)
    df <- left_join(df, cordat[, c("SampleID", "gdm", "birth_weight", "sex", "GA", "mother_ethnicity")], by = "SampleID")
    df$mother_ethnicity <- as.numeric(as.factor(df$mother_ethnicity))
    
    # Model components
    direct_effect <- "
  # Direct effect
  birth_weight ~ c*gdm + sex_bw*sex + GA_bw*GA + ethnicity_bw*mother_ethnicity
"
    mediation_model <- paste0(direct_effect, "\n", generate_mediation_model_body(i))
    
    # Perform mediation analysis
    fit <- sem(mediation_model, data = df)
    param_est <- parameterEstimates(fit, standardized = TRUE, ci = TRUE)
    
    # Extract total_indirect statistics
    ind_row <- param_est %>% filter(label == "total_indirect") %>%
      select(est = est, ci.lower = ci.lower, ci.upper = ci.upper, pvalue = pvalue)
    
    # Save results
    results[[i]] <- data.frame(
      PCs = i,
      percent_variance = round(pca_cumvar[i] * 100, 2),
      total_indirect = ind_row$est,
      CI_lower = ind_row$ci.lower,
      CI_upper = ind_row$ci.upper,
      p = ind_row$pvalue
    )
  }
  
  final_results <- bind_rows(results)
  return(final_results)
}

################################################################################
# GUSTO Principal Component Mediation Analysis
################################################################################

## GUSTO assembly: PC-based mediation analysis
GUSTO.expr.assembly <- data.frame(GUSTO.vst.assembly)
GUSTO.expr.assembly$transcript_id <- row.names(GUSTO.expr.assembly)

# Filter to most significant transcripts (p < 0.00005) for PC analysis
GUSTO.expr.assembly <- GUSTO.expr.assembly[GUSTO.expr.assembly$transcript_id %in% 
                                             results.assembly.GUSTO.gdm[results.assembly.GUSTO.gdm$pvalue < 0.00005, "transcript_id"], ]
GUSTO.expr.assembly$transcript_id <- NULL

## Examine principal components 1-10
pca_result <- prcomp(t(as.matrix(GUSTO.expr.assembly)))
variance_explained <- pca_result$sdev^2
percent_variance <- variance_explained / sum(variance_explained) * 100
percent_variance[1:10]
cumulative_variance <- cumsum(percent_variance)
cumulative_variance[1:10]
cat("PC1-5 individual variance explained:\n")
round(percent_variance[1:5], 2)
cat("\nPC1-5 cumulative variance explained:\n")
round(cumulative_variance[1:5], 2)
cat("\nTotal variance explained by PC1-5:", round(cumulative_variance[5], 2), "%\n")

plot(percent_variance[1:20], type = "b", 
     xlab = "Principal Component", 
     ylab = "Percentage Variance Explained")

## Examine correlation with maternal ethnicity (potential confounder)
pca_result <- prcomp(t(as.matrix(GUSTO.expr.assembly)))
pc_scores <- data.frame(pca_result$x[, 1:10])
pc_scores$SampleID <- rownames(pc_scores)
pc_ethnicity <- merge(pc_scores, cordat.GUSTO[, c("SampleID", "mother_ethnicity")], by = "SampleID")

for(i in 1:10) {
  pc_name <- paste0("PC", i)
  formula_str <- paste(pc_name, "~ mother_ethnicity")
  aov_result <- aov(as.formula(formula_str), data = pc_ethnicity)
  p_value <- summary(aov_result)[[1]][["Pr(>F)"]][1]
  
  # Suppress the message about eta squared
  eta_sq <- suppressMessages(eta_squared(aov_result)$Eta2)
  
  cat(pc_name, "vs ethnicity p-value:", p_value, 
      ", eta-squared:", round(eta_sq, 4), "\n")
}

# Run PC sensitivity analysis for GUSTO assembly
GUSTO.assembly.PCsensitivity <- PCmed.sensitivity(GUSTO.expr.assembly, cordat.GUSTO, 10)

## GUSTO gencode: PC-based mediation analysis
GUSTO.expr.gencode <- data.frame(GUSTO.vst.gencode)
GUSTO.expr.gencode$transcript_id <- row.names(GUSTO.expr.gencode)

# Filter to most significant transcripts for PC analysis
GUSTO.expr.gencode <- GUSTO.expr.gencode[GUSTO.expr.gencode$transcript_id %in% 
                                           results.gencode.GUSTO.gdm[results.gencode.GUSTO.gdm$pvalue < 0.00005, "transcript_id"], ]
GUSTO.expr.gencode$transcript_id <- NULL

## Examine principal components 1-10
pca_result <- prcomp(t(as.matrix(GUSTO.expr.gencode)))
variance_explained <- pca_result$sdev^2
percent_variance <- variance_explained / sum(variance_explained) * 100
percent_variance[1:10]
cumulative_variance <- cumsum(percent_variance)
cumulative_variance[1:10] # 8 PCs = cumulative variance of 76%
cat("PC1-8 individual variance explained:\n")
round(percent_variance[1:8], 2)
cat("\nPC1-8 cumulative variance explained:\n")
round(cumulative_variance[1:8], 2)
cat("\nTotal variance explained by PC1-8:", round(cumulative_variance[8], 2), "%\n")

plot(percent_variance[1:20], type = "b", 
     xlab = "Principal Component", 
     ylab = "Percentage Variance Explained")

## Examine correlation with maternal ethnicity
pca_result <- prcomp(t(as.matrix(GUSTO.expr.gencode)))
pc_scores <- data.frame(pca_result$x[, 1:10])
pc_scores$SampleID <- rownames(pc_scores)
pc_ethnicity <- merge(pc_scores, cordat.GUSTO[, c("SampleID", "mother_ethnicity")], by = "SampleID")

for(i in 1:10) {
  pc_name <- paste0("PC", i)
  formula_str <- paste(pc_name, "~ mother_ethnicity")
  aov_result <- aov(as.formula(formula_str), data = pc_ethnicity)
  p_value <- summary(aov_result)[[1]][["Pr(>F)"]][1]
  
  # Suppress the message about eta squared
  eta_sq <- suppressMessages(eta_squared(aov_result)$Eta2)
  
  cat(pc_name, "vs ethnicity p-value:", p_value, 
      ", eta-squared:", round(eta_sq, 4), "\n")
}

# Run PC sensitivity analysis for GUSTO gencode
GUSTO.gencode.PCsensitivity <- PCmed.sensitivity(GUSTO.expr.gencode, cordat.GUSTO, 10)

## Create sensitivity plot for GUSTO (Figure 4B)
GUSTO.assembly.PCsensitivity$Assembly <- "LR"
GUSTO.gencode.PCsensitivity$Assembly <- "GENCODEv45"
PCsensitivity <- rbind(GUSTO.assembly.PCsensitivity,
                       GUSTO.gencode.PCsensitivity)
rownames(PCsensitivity) <- NULL

# Figure 4B - PC mediation sensitivity analysis for GUSTO
fig4b <- ggplot(PCsensitivity, aes(x = PCs, color = Assembly, fill = Assembly)) +
  geom_hline(yintercept = 0, linetype = "dotted") + 
  geom_point(aes(y = total_indirect), size = 1.5, position = position_dodge(width = 0.5)) +
  geom_line(aes(y = total_indirect, group = Assembly, color = Assembly), position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = CI_lower, ymax = CI_upper, color = Assembly), position = position_dodge(width = 0.5)) +
  # Add percent_variance on secondary y-axis
  geom_point(aes(y = percent_variance/100 * diff(range(PCsensitivity$total_indirect)) + min(PCsensitivity$total_indirect)),
             shape = 17, size = 2, alpha = 0.7, position = position_dodge(width = 0.5)) +
  geom_line(aes(y = percent_variance/100 * diff(range(PCsensitivity$total_indirect)) + min(PCsensitivity$total_indirect),
                group = Assembly), linetype = "dashed", alpha = 0.7, position = position_dodge(width = 0.5)) +
  scale_x_continuous(breaks = seq.int(1, 10, 1)) +
  # Add secondary y-axis
  scale_y_continuous(
    name = "Indirect Effect Estimate",
    sec.axis = sec_axis(~ (. - min(PCsensitivity$total_indirect)) / diff(range(PCsensitivity$total_indirect)) * 100,
                        name = "% Variance Explained")
  ) +
  labs(title = "Assembly Estimates with Confidence Bands and % Variance Explained",
       x = "PCs (Principal Components)",
       color = "Assembly",
       fill = "Assembly") +
  theme_minimal() +
  theme(legend.position = "bottom", title = element_text(size = 8)) +
  scale_fill_manual(values = c("royalblue", "seagreen")) +
  scale_color_manual(values = c("royalblue", "seagreen"))

fig4b
ggsave("fig4b.png", fig4b, width = 5.10, height = 3.45, units = "in", dpi = 300)

################################################################################
# Gen3G Principal Component Mediation Analysis
################################################################################

# Function for Gen3G PC sensitivity (simpler model without ethnicity)
PCmed.sensitivity.Gen3G <- function(expmat, cordat, PCmax) {
  
  require(lavaan)
  require(dplyr)
  
  generate_mediation_model_body <- function(n_pcs) {
    mediator_block <- paste0(
      "  PC", 1:n_pcs, " ~ a", 1:n_pcs, "*gdm + sex_m", 1:n_pcs, "*sex + GA_m", 1:n_pcs, "*GA",
      collapse = "\n"
    )
    
    outcome_block <- paste0(
      "  birth_weight ~ ", paste0("b", 1:n_pcs, "*PC", 1:n_pcs, collapse = " + ")
    )
    
    indirect_block <- paste0(
      "  ind", 1:n_pcs, " := a", 1:n_pcs, "*b", 1:n_pcs, collapse = "\n"
    )
    
    total_indirect_block <- paste0(
      "  total_indirect := ", paste0("ind", 1:n_pcs, collapse = " + ")
    )
    
    total_effect_block <- "  total := c + total_indirect"
    prop_mediated_block <- "  prop_mediated := total_indirect/total"
    
    prop_each_block <- paste0(
      "  prop_med", 1:n_pcs, " := ind", 1:n_pcs, "/total", collapse = "\n"
    )
    
    paste(
      "# Mediator paths with sex and GA adjustment\n", mediator_block,
      "\n\n# Outcome paths\n", outcome_block,
      "\n\n# Compute indirect effects\n", indirect_block,
      "\n\n# Total indirect effect\n", total_indirect_block,
      "\n\n# Total effect\n", total_effect_block,
      "\n\n# Proportion mediated (overall)\n", prop_mediated_block,
      "\n\n# Proportion for each mediator\n", prop_each_block,
      "\n", sep = ""
    )
  }
  
  pca_res <- prcomp(t(as.matrix(expmat)))
  pca_df <- as.data.frame(pca_res$x[, 1:PCmax])
  pca_var_explained <- summary(pca_res)$importance["Proportion of Variance", 1:PCmax]
  pca_cumvar <- cumsum(pca_var_explained)
  
  results <- list()
  
  for (i in 1:PCmax) {
    df <- data.frame(pca_df[, 1:i])
    colnames(df) <- paste0("PC", 1:i)
    df$Run <- rownames(pca_df)
    df <- left_join(df, cordat[, c("Run", "gdm", "birth_weight", "sex", "GA")], by = "Run")
    
    direct_effect <- "
  # Direct effect
  birth_weight ~ c*gdm + sex_bw*sex + GA_bw*GA
"
    mediation_model <- paste0(direct_effect, "\n", generate_mediation_model_body(i))
    
    fit <- sem(mediation_model, data = df)
    param_est <- parameterEstimates(fit, standardized = TRUE, ci = TRUE)
    
    ind_row <- param_est %>% filter(label == "total_indirect") %>%
      select(est = est, ci.lower = ci.lower, ci.upper = ci.upper, pvalue = pvalue)
    
    results[[i]] <- data.frame(
      PCs = i,
      percent_variance = round(pca_cumvar[i] * 100, 2),
      total_indirect = ind_row$est,
      CI_lower = ind_row$ci.lower,
      CI_upper = ind_row$ci.upper,
      p = ind_row$pvalue
    )
  }
  
  final_results <- bind_rows(results)
  return(final_results)
}

## Gen3G assembly: PC-based mediation analysis
Gen3G.expr.assembly <- data.frame(Gen3G.vst.assembly)
Gen3G.expr.assembly$transcript_id <- row.names(Gen3G.expr.assembly)

# Filter to most significant transcripts for PC analysis
Gen3G.expr.assembly <- Gen3G.expr.assembly[Gen3G.expr.assembly$transcript_id %in% 
                                             results.assembly.Gen3G.gdm[results.assembly.Gen3G.gdm$pvalue < 0.00005, "transcript_id"], ]
Gen3G.expr.assembly$transcript_id <- NULL

## Examine principal components 1-20
pca_result <- prcomp(t(as.matrix(Gen3G.expr.assembly)))
variance_explained <- pca_result$sdev^2
percent_variance <- variance_explained / sum(variance_explained) * 100
percent_variance[1:20]
cumulative_variance <- cumsum(percent_variance)
cumulative_variance[1:20] # 7 PCs = 77.39% variance

plot(percent_variance[1:20], type = "b", 
     xlab = "Principal Component", 
     ylab = "Percentage Variance Explained")

# Run PC sensitivity analysis for Gen3G assembly
Gen3G.assembly.PCsensitivity <- PCmed.sensitivity.Gen3G(Gen3G.expr.assembly, cordat.Gen3G, 20)

## Gen3G gencode: PC-based mediation analysis
Gen3G.expr.gencode <- data.frame(Gen3G.vst.gencode)
Gen3G.expr.gencode$transcript_id <- row.names(Gen3G.expr.gencode)

# Filter to most significant transcripts for PC analysis
Gen3G.expr.gencode <- Gen3G.expr.gencode[Gen3G.expr.gencode$transcript_id %in% 
                                           results.gencode.Gen3G.gdm[results.gencode.Gen3G.gdm$pvalue < 0.00005, "transcript_id"], ]
Gen3G.expr.gencode$transcript_id <- NULL

## Examine principal components 1-20
pca_result <- prcomp(t(as.matrix(Gen3G.expr.gencode)))
variance_explained <- pca_result$sdev^2
percent_variance <- variance_explained / sum(variance_explained) * 100
percent_variance[1:20]
cumulative_variance <- cumsum(percent_variance)
cumulative_variance[1:20] # 14 PCs = cumulative variance of 75.9%

plot(percent_variance[1:20], type = "b", 
     xlab = "Principal Component", 
     ylab = "Percentage Variance Explained")

# Run PC sensitivity analysis for Gen3G gencode
Gen3G.gencode.PCsensitivity <- PCmed.sensitivity.Gen3G(Gen3G.expr.gencode, cordat.Gen3G, 20)

## Create sensitivity plot for Gen3G (Figure S8)
Gen3G.assembly.PCsensitivity$Assembly <- "LR"
Gen3G.gencode.PCsensitivity$Assembly <- "GENCODEv45"
PCsensitivity <- rbind(Gen3G.assembly.PCsensitivity,
                       Gen3G.gencode.PCsensitivity)
rownames(PCsensitivity) <- NULL

# Figure S8 - PC mediation sensitivity analysis for Gen3G
figs8 <- ggplot(PCsensitivity, aes(x = PCs, color = Assembly, fill = Assembly)) +
  geom_hline(yintercept = 0, linetype = "dotted") + 
  geom_point(aes(y = total_indirect), size = 2, position = position_dodge(width = 0.5)) +
  geom_line(aes(y = total_indirect, group = Assembly, color = Assembly), position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = CI_lower, ymax = CI_upper, color = Assembly), position = position_dodge(width = 0.5)) +
  geom_point(aes(y = percent_variance/100 * diff(range(PCsensitivity$total_indirect)) + min(PCsensitivity$total_indirect)),
             shape = 17, size = 2.5, alpha = 0.7, position = position_dodge(width = 0.5)) +
  geom_line(aes(y = percent_variance/100 * diff(range(PCsensitivity$total_indirect)) + min(PCsensitivity$total_indirect),
                group = Assembly), linetype = "dashed", alpha = 0.7, position = position_dodge(width = 0.5)) +
  scale_x_continuous(breaks = seq.int(1, 20, 1)) +
  scale_y_continuous(
    name = "Indirect Effect Estimate",
    sec.axis = sec_axis(~ (. - min(PCsensitivity$total_indirect)) / diff(range(PCsensitivity$total_indirect)) * 100,
                        name = "% Variance Explained")
  ) +
  labs(title = "Assembly Estimates with Confidence Bands and % Variance Explained",
       x = "PCs (Principal Components)",
       color = "Assembly",
       fill = "Assembly") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("royalblue", "seagreen")) +
  scale_color_manual(values = c("royalblue", "seagreen"))

figs8

################################################################################
# Gene Ontology Enrichment Analysis Setup
################################################################################

# Create directories for enrichment analysis results
dir.create("fig5_SC/EnrichR_Results_DTE", showWarnings = FALSE)
dir.create("fig5_SC/EnrichR_Results_DTE/Plots", showWarnings = FALSE)
dir.create("fig5_SC/EnrichR_Results_DGE", showWarnings = FALSE)
dir.create("fig5_SC/EnrichR_Results_DGE/Plots", showWarnings = FALSE)

# Function to convert Ensembl IDs to gene symbols for enrichR compatibility
ensembl_to_symbol <- function(ensembl_ids) {
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

# Function to run enrichR analysis and save all results
run_enrichr_analysis <- function(gene_list, background_list, name, databases, dir_out) {
  
  if (length(gene_list) < 5) {
    cat("Skipping", name, "analysis: too few genes\n")
    return(NULL)
  }
  
  all_results <- list()
  
  tryCatch({
    # Run enrichR with background list
    enrichr_results <- enrichr(
      genes = gene_list,
      databases = databases,
      background = background_list
    )
    
    # Process each database result
    for (db in names(enrichr_results)) {
      result_df <- enrichr_results[[db]]
      
      # Skip if no results
      if (is.null(result_df) || nrow(result_df) == 0) {
        cat("No results found for", name, "in database:", db, "\n")
        next
      }
      
      # Add dataset info
      result_df$Dataset <- name
      
      # Save all results regardless of p-value
      write.csv(
        result_df,
        file = paste0(dir_out, "/", name, "_", gsub("[^a-zA-Z0-9]", "_", db), ".csv"),
        row.names = FALSE)
      
      cat("Saved", nrow(result_df), "terms for", name, "in database:", db, "\n")
      
      # Store in all_results
      all_results[[db]] <- result_df
    }
    
    return(all_results)
  }, error = function(e) {
    cat("Error in enrichR analysis for", name, ":", conditionMessage(e), "\n")
    return(NULL)
  })
}

# Create combined summary table function
create_summary_table <- function(enrichr_assembly_results,
                                 enrichr_gencode_results,
                                 file_out) {
  # Initialize empty data frame
  all_results <- data.frame()
  
  # Process LR-assembly results
  for (db in dbs) {
    if (!is.null(enrichr_assembly_results[[db]])) {
      df <- enrichr_assembly_results[[db]]
      df$Database <- db
      df$Gene_Set <- "LR-assembly_significant"
      df$Background_Set <- "LR-assembly_universe"
      all_results <- rbind(all_results, df)
    }
  }
  
  # Process GENCODE results
  for (db in dbs) {
    if (!is.null(enrichr_gencode_results[[db]])) {
      df <- enrichr_gencode_results[[db]]
      df$Database <- db
      df$Gene_Set <- "GENCODE_significant"
      df$Background_Set <- "GENCODE_universe"
      all_results <- rbind(all_results, df)
    }
  }
  
  if (nrow(all_results) > 0) {
    write.csv(all_results, file_out, row.names = FALSE)
    cat("Created summary table with", nrow(all_results), "terms\n")
  } else {
    cat("No terms found to create summary table\n")
  }
  
  return(all_results)
}

# Function to check if a plot could benefit from an axis break
needs_axis_break <- function(y, threshold_ratio = 2, top_fraction = 0.1, buffer = 1) {
  y <- y[!is.na(y)]
  if (length(y) < 5) return(c(FALSE, NA, NA))
  
  y_sorted <- sort(y, decreasing = TRUE)
  n_top <- max(1, floor(length(y) * top_fraction))
  
  top_values <- y_sorted[1:n_top]
  rest_values <- y_sorted[-(1:n_top)]
  
  if (length(rest_values) == 0) return(c(FALSE, NA, NA))
  
  ratio <- mean(top_values) / mean(rest_values)
  
  if (ratio >= threshold_ratio) {
    lower_break <- max(rest_values) + buffer
    upper_break <- min(top_values) - buffer
    if (lower_break >= upper_break) {
      lower_break <- max(rest_values) + 0.5
      upper_break <- min(top_values) - 0.5
    }
    return(c(TRUE, lower_break, upper_break))
  } else {
    return(c(FALSE, NA, NA))
  }
}

# Function to create lollipop plot showing both p-value and odds ratio
create_lollipop_plot <- function(data, dataset_name, n = 15, output_file) {
  # Filter to specific dataset and sort by adjusted p-value
  df <- data %>% 
    filter(Dataset == dataset_name) %>%
    arrange(Adjusted.P.value) %>%
    head(n)
  
  # Create abbreviated ontology labels for cleaner display
  df$Ontology <- df$Term
  df$Ontology <- ifelse(
    nchar(df$Ontology) > 60,
    paste0(substr(df$Ontology, 1, 57), "..."),
    df$Ontology
  )
  
  # Calculate -log10(p-value) for visualization
  df$neglog10p <- -log10(df$Adjusted.P.value)
  
  # Reorder factors for the plot
  df$Ontology <- fct_reorder(df$Ontology, df$neglog10p)
  
  df <- df %>%
    mutate(Database = recode(Database,
                             "GO_Biological_Process_2023" = "GOBP",
                             "GO_Molecular_Function_2023" = "GOMF",
                             "GO_Cellular_Component_2023" = "GOCC",
                             "KEGG_2021_Human" = "KEGG",
                             "WikiPathway_2021_Human" = "WikiPathway",
                             "Reactome_2022" = "Reactome",
                             "MSigDB_Hallmark_2020" = "MSigDB_Hallmark"
    ))
  
  # Define color palette based on database
  db_palette <- c(
    "GOBP" = "#E41A1C",
    "GOMF" = "#377EB8",
    "GOCC" = "#4DAF4A",
    "KEGG" = "#984EA3",
    "WikiPathway" = "#FF7F00",
    "Reactome" = "#FFFF33",
    "MSigDB_Hallmark" = "#A65628"
  )
  
  # Generate missing colors for any databases not in our palette
  all_dbs <- unique(df$Database)
  missing_dbs <- setdiff(all_dbs, names(db_palette))
  
  if (length(missing_dbs) > 0) {
    new_colors <- rainbow(length(missing_dbs))
    names(new_colors) <- missing_dbs
    db_palette <- c(db_palette, new_colors)
  }
  
  # Create the lollipop plot
  p <- ggplot(df, aes(x = neglog10p, y = Ontology)) +
    # Add segments for the "stick"
    geom_segment(
      aes(x = 0, xend = neglog10p, y = Ontology, yend = Ontology),
      color = "gray50", linewidth = 0.7
    ) +
    # Add points for the "lollipop" - sized by odds ratio, colored by database
    geom_point(
      aes(size = Odds.Ratio, color = Database),
      alpha = 0.8
    ) +
    # Scale settings
    scale_size_continuous(range = c(3, 8), name = "Odds Ratio") +
    scale_color_manual(values = db_palette) +
    # Labels and theme
    labs(
      title = paste("Top 15 Enriched Terms -", dataset_name),
      subtitle = "Point size = Odds Ratio, X-axis = -log10(FDR)",
      x = "-log10(Adjusted P-value)",
      y = ""
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10),
      axis.text.y = element_text(size = 8),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.title = element_blank()
    )
  
  # Add text labels for odds ratio
  if ("Overlap" %in% colnames(df)) {
    p <- p + geom_text(
      aes(label = sprintf("OR=%.1f (%s)", Odds.Ratio, Overlap)),
      hjust = -0.2, vjust = 0.5, size = 2.5
    )
  } else {
    p <- p + geom_text(
      aes(label = sprintf("OR=%.1f", Odds.Ratio)),
      hjust = -0.2, vjust = 0.5, size = 2.5
    )
  }
  
  # Check if the plot needs an axis break
  cab <- needs_axis_break(df$neglog10p)
  if(cab[1] == 1 & ((cab[3] - cab[2]) > 3)) {
    p <- p + scale_x_break(c(cab[2], cab[3]), expand = F, space = 0) +
      xlim(0, max(df$neglog10p + 1)) +
      theme(
        axis.ticks.x.top = element_blank(),
        axis.text.x.top = element_blank())
  } else {
    p <- p +
      scale_x_continuous(expand = expansion(mult = 0),
                         limits = c(0, max(df$neglog10p + 1)))
  }
  
  # Save the plot
  invisible(ggsave(output_file, p, width = 6, height = 8))
  
  return(p)
}

# Function for creating lollipop plots for each dataset within a results table
make_lolipops_perdataset <- function(datasets, all_results, dir_out, file_suffix) {
  for (dataset in datasets) {
    # Clean dataset name for file naming
    clean_name <- gsub("[^a-zA-Z0-9]", "_", dataset)
    
    # Create and save the lollipop plot
    p <- create_lollipop_plot(
      all_results, 
      dataset, 
      15, 
      paste0(dir_out, clean_name, file_suffix)
    )
  }
}

# Set up enrichR databases
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

################################################################################
# Gene Universe Preparation
################################################################################

# Extract gene ID from tx2g tables (remove version number)
tx2g.gencode$gene_id <- sapply(strsplit(tx2g.gencode$gene_id, '[.]'),
                               function(x) x[1])
tx2g.assembly$gene_id <- sapply(strsplit(tx2g.assembly$gene_id, '[.]'),
                                function(x) x[1])

# Define universes for each annotation
## For GENCODE: all genes
universe_gencode <- unique(tx2g.gencode$gene_id)
## For LR-assembly: only genes that start with ENSG
universe_assembly <- unique(tx2g.assembly$gene_id[grepl("^ENSG", tx2g.assembly$gene_id)])

# Convert universes to symbols
universe_assembly_result <- ensembl_to_symbol(universe_assembly)
symbols_universe_assembly <- universe_assembly_result$symbols
mapping_universe_assembly <- universe_assembly_result$mapping

universe_gencode_result <- ensembl_to_symbol(universe_gencode)
symbols_universe_gencode <- universe_gencode_result$symbols
mapping_universe_gencode <- universe_gencode_result$mapping

################################################################################
# Top 50 Isoform-Rich Genes Analysis
################################################################################

# Identify 50 most isoform-rich genes for pathway analysis
iso.tab <- data.frame(table(classif$associated_gene))
iso.tab <- iso.tab[rev(order(iso.tab$Freq)), ]
top50 <- iso.tab[1:50, 1]
top50 <- sapply(strsplit(as.character(top50), '[.]'), function(x) x[1])

top50 <- ensembl_to_symbol(top50)
symbols_top50 <- top50$symbols
mapping_top50 <- top50$mapping

# Run enrichment analysis on top 50 isoform-rich genes
top50_enrichr <- run_enrichr_analysis(symbols_top50, 
                                      symbols_universe_assembly, 
                                      "top50isorichgenes", dbs,
                                      "fig5_SC/EnrichR_Results_top50isorichgenes")
top50_summary <- create_summary_table(top50_enrichr,
                                      NULL,
                                      "fig5_SC/EnrichR_Results_top50isorichgenes/results.csv")

# Extract molecular function results for top 50 isoform-rich genes
top50MF <- top50_summary[top50_summary$Database == "GO_Molecular_Function_2023", ]
write.csv(top50MF, "top50isorichgenes_GOMF.csv", row.names = FALSE)

################################################################################
# Differential Transcript Expression (DTE) Enrichment Analysis
################################################################################

## Prepare transcript results and map to genes
# Set first column to transcript_id in each dataset
results.assembly.Gen3G.gdm <- results.assembly.Gen3G.gdm[, c(7, 1:6)]
results.gencode.Gen3G.gdm <- results.gencode.Gen3G.gdm[, c(7, 1:6)]
results.assembly.GUSTO.gdm <- results.assembly.GUSTO.gdm[, c(7, 1:6)]
results.gencode.GUSTO.gdm <- results.gencode.GUSTO.gdm[, c(7, 1:6)]

## Identify transcripts with padj < 0.1
sig_transcripts_assembly_Gen3G <- results.assembly.Gen3G.gdm[results.assembly.Gen3G.gdm$padj < 0.1, ]
sig_transcripts_gencode_Gen3G <- results.gencode.Gen3G.gdm[results.gencode.Gen3G.gdm$padj < 0.1, ]
sig_transcripts_assembly_GUSTO <- results.assembly.GUSTO.gdm[results.assembly.GUSTO.gdm$padj < 0.1, ]
sig_transcripts_gencode_GUSTO <- results.gencode.GUSTO.gdm[results.gencode.GUSTO.gdm$padj < 0.1, ]

## Map transcripts to genes
sig_transcripts_assembly_Gen3G <- left_join(sig_transcripts_assembly_Gen3G, tx2g.assembly)
sig_transcripts_gencode_Gen3G <- left_join(sig_transcripts_gencode_Gen3G, tx2g.gencode)
sig_transcripts_assembly_GUSTO <- left_join(sig_transcripts_assembly_GUSTO, tx2g.assembly)
sig_transcripts_gencode_GUSTO <- left_join(sig_transcripts_gencode_GUSTO, tx2g.gencode)

## Get unique gene IDs
unique_sig_genes_assembly_Gen3G <- unique(sig_transcripts_assembly_Gen3G$gene_id)
unique_sig_genes_gencode_Gen3G <- unique(sig_transcripts_gencode_Gen3G$gene_id)
unique_sig_genes_assembly_GUSTO <- unique(sig_transcripts_assembly_GUSTO$gene_id)
unique_sig_genes_gencode_GUSTO <- unique(sig_transcripts_gencode_GUSTO$gene_id)

## Filter significant assembly genes to only include ENSG genes
unique_sig_genes_assembly_Gen3G_filtered <- unique_sig_genes_assembly_Gen3G[grepl("^ENSG", unique_sig_genes_assembly_Gen3G)]
unique_sig_genes_assembly_GUSTO_filtered <- unique_sig_genes_assembly_GUSTO[grepl("^ENSG", unique_sig_genes_assembly_GUSTO)]

## Convert gene sets to symbols
assembly_result_Gen3G <- ensembl_to_symbol(unique_sig_genes_assembly_Gen3G_filtered)
symbols_sig_assembly_Gen3G <- assembly_result_Gen3G$symbols
mapping_sig_assembly_Gen3G <- assembly_result_Gen3G$mapping

gencode_result_Gen3G <- ensembl_to_symbol(unique_sig_genes_gencode_Gen3G)
symbols_sig_gencode_Gen3G <- gencode_result_Gen3G$symbols
mapping_sig_gencode_Gen3G <- gencode_result_Gen3G$mapping

assembly_result_GUSTO <- ensembl_to_symbol(unique_sig_genes_assembly_GUSTO_filtered)
symbols_sig_assembly_GUSTO <- assembly_result_GUSTO$symbols
mapping_sig_assembly_GUSTO <- assembly_result_GUSTO$mapping

gencode_result_GUSTO <- ensembl_to_symbol(unique_sig_genes_gencode_GUSTO)
symbols_sig_gencode_GUSTO <- gencode_result_GUSTO$symbols
mapping_sig_gencode_GUSTO <- gencode_result_GUSTO$mapping

## Run enrichR analysis for DTE
enrichr_assembly_results_Gen3G <- run_enrichr_analysis(symbols_sig_assembly_Gen3G, 
                                                       symbols_universe_assembly, 
                                                       "Gen3G_assembly_sig", dbs,
                                                       "fig5_SC/EnrichR_Results_DTE")

enrichr_gencode_results_Gen3G <- run_enrichr_analysis(symbols_sig_gencode_Gen3G, 
                                                      symbols_universe_gencode, 
                                                      "Gen3G_gencode_sig", dbs,
                                                      "fig5_SC/EnrichR_Results_DTE")

enrichr_assembly_results_GUSTO <- run_enrichr_analysis(symbols_sig_assembly_GUSTO, 
                                                       symbols_universe_assembly, 
                                                       "GUSTO_assembly_sig", dbs,
                                                       "fig5_SC/EnrichR_Results_DTE")

enrichr_gencode_results_GUSTO <- run_enrichr_analysis(symbols_sig_gencode_GUSTO, 
                                                      symbols_universe_gencode, 
                                                      "GUSTO_gencode_sig", dbs,
                                                      "fig5_SC/EnrichR_Results_DTE")

## Create combined summary tables for DTE
all_results_Gen3G <- create_summary_table(enrichr_assembly_results_Gen3G,
                                          enrichr_gencode_results_Gen3G,
                                          "fig5_SC/EnrichR_Results_DTE/results_Gen3G_DTE.csv")
all_results_Gen3G.DTE <- all_results_Gen3G

all_results_GUSTO <- create_summary_table(enrichr_assembly_results_GUSTO,
                                          enrichr_gencode_results_GUSTO,
                                          "fig5_SC/EnrichR_Results_DTE/results_GUSTO_DTE.csv")
all_results_GUSTO.DTE <- all_results_GUSTO

## Create lollipop plots for DTE visualization
# Get list of datasets
datasets_Gen3G <- unique(all_results_Gen3G$Dataset)
datasets_GUSTO <- unique(all_results_GUSTO$Dataset)

# Create lollipop plots for each dataset
make_lolipops_perdataset(datasets_Gen3G, all_results_Gen3G,
                         "fig5_SC/EnrichR_Results_DTE/Plots/",
                         "_lollipop_DTE.pdf")

make_lolipops_perdataset(datasets_GUSTO, all_results_GUSTO,
                         "fig5_SC/EnrichR_Results_DTE/Plots/",
                         "_lollipop_DTE.pdf")

################################################################################
# Differential Gene Expression (DGE) Enrichment Analysis
################################################################################

## Prepare gene-level results
# Subset columns and set Confirmation.P as padj
results.assembly.Gen3G.gdm.g <- results.assembly.Gen3G.gdm.g[, c(1:9)]
names(results.assembly.Gen3G.gdm.g)[9] <- "padj"
results.assembly.Gen3G.gdm.g$gene_id <- sapply(strsplit(results.assembly.Gen3G.gdm.g$gene_id, '[.]'), function(x) x[1])

results.gencode.Gen3G.gdm.g <- results.gencode.Gen3G.gdm.g[, c(1:9)]
names(results.gencode.Gen3G.gdm.g)[9] <- "padj"
results.gencode.Gen3G.gdm.g$gene_id <- sapply(strsplit(results.gencode.Gen3G.gdm.g$gene_id, '[.]'), function(x) x[1])

results.assembly.GUSTO.gdm.g <- results.assembly.GUSTO.gdm.g[, c(1:9)]
names(results.assembly.GUSTO.gdm.g)[9] <- "padj"
results.assembly.GUSTO.gdm.g$gene_id <- sapply(strsplit(results.assembly.GUSTO.gdm.g$gene_id, '[.]'), function(x) x[1])

results.gencode.GUSTO.gdm.g <- results.gencode.GUSTO.gdm.g[, c(1:9)]
names(results.gencode.GUSTO.gdm.g)[9] <- "padj"
results.gencode.GUSTO.gdm.g$gene_id <- sapply(strsplit(results.gencode.GUSTO.gdm.g$gene_id, '[.]'), function(x) x[1])

## Identify genes with padj < 0.1
sig_transcripts_assembly_Gen3G <- results.assembly.Gen3G.gdm.g[results.assembly.Gen3G.gdm.g$padj < 0.1, ]
sig_transcripts_gencode_Gen3G <- results.gencode.Gen3G.gdm.g[results.gencode.Gen3G.gdm.g$padj < 0.1, ]
sig_transcripts_assembly_GUSTO <- results.assembly.GUSTO.gdm.g[results.assembly.GUSTO.gdm.g$padj < 0.1, ]
sig_transcripts_gencode_GUSTO <- results.gencode.GUSTO.gdm.g[results.gencode.GUSTO.gdm.g$padj < 0.1, ]

## Get unique gene IDs
unique_sig_genes_assembly_Gen3G <- unique(sig_transcripts_assembly_Gen3G$gene_id)
unique_sig_genes_gencode_Gen3G <- unique(sig_transcripts_gencode_Gen3G$gene_id)
unique_sig_genes_assembly_GUSTO <- unique(sig_transcripts_assembly_GUSTO$gene_id)
unique_sig_genes_gencode_GUSTO <- unique(sig_transcripts_gencode_GUSTO$gene_id)

## Filter significant assembly genes to only include ENSG genes
unique_sig_genes_assembly_Gen3G_filtered <- unique_sig_genes_assembly_Gen3G[grepl("^ENSG", unique_sig_genes_assembly_Gen3G)]
unique_sig_genes_assembly_GUSTO_filtered <- unique_sig_genes_assembly_GUSTO[grepl("^ENSG", unique_sig_genes_assembly_GUSTO)]

## Convert gene sets to symbols
assembly_result_Gen3G <- ensembl_to_symbol(unique_sig_genes_assembly_Gen3G_filtered)
symbols_sig_assembly_Gen3G <- assembly_result_Gen3G$symbols
mapping_sig_assembly_Gen3G <- assembly_result_Gen3G$mapping

gencode_result_Gen3G <- ensembl_to_symbol(unique_sig_genes_gencode_Gen3G)
symbols_sig_gencode_Gen3G <- gencode_result_Gen3G$symbols
mapping_sig_gencode_Gen3G <- gencode_result_Gen3G$mapping

assembly_result_GUSTO <- ensembl_to_symbol(unique_sig_genes_assembly_GUSTO_filtered)
symbols_sig_assembly_GUSTO <- assembly_result_GUSTO$symbols
mapping_sig_assembly_GUSTO <- assembly_result_GUSTO$mapping

gencode_result_GUSTO <- ensembl_to_symbol(unique_sig_genes_gencode_GUSTO)
symbols_sig_gencode_GUSTO <- gencode_result_GUSTO$symbols
mapping_sig_gencode_GUSTO <- gencode_result_GUSTO$mapping

## Run enrichR analysis for DGE
enrichr_assembly_results_Gen3G <- run_enrichr_analysis(symbols_sig_assembly_Gen3G, 
                                                       symbols_universe_assembly, 
                                                       "Gen3G_assembly_sig", dbs,
                                                       "fig5_SC/EnrichR_Results_DGE")

enrichr_gencode_results_Gen3G <- run_enrichr_analysis(symbols_sig_gencode_Gen3G, 
                                                      symbols_universe_gencode, 
                                                      "Gen3G_gencode_sig", dbs,
                                                      "fig5_SC/EnrichR_Results_DGE")

enrichr_assembly_results_GUSTO <- run_enrichr_analysis(symbols_sig_assembly_GUSTO, 
                                                       symbols_universe_assembly, 
                                                       "GUSTO_assembly_sig", dbs,
                                                       "fig5_SC/EnrichR_Results_DGE")

enrichr_gencode_results_GUSTO <- run_enrichr_analysis(symbols_sig_gencode_GUSTO, 
                                                      symbols_universe_gencode, 
                                                      "GUSTO_gencode_sig", dbs,
                                                      "fig5_SC/EnrichR_Results_DGE")

## Create combined summary tables for DGE
all_results_Gen3G <- create_summary_table(enrichr_assembly_results_Gen3G,
                                          enrichr_gencode_results_Gen3G,
                                          "fig5_SC/EnrichR_Results_DGE/results_Gen3G_DGE.csv")
all_results_Gen3G.DGE <- all_results_Gen3G

all_results_GUSTO <- create_summary_table(enrichr_assembly_results_GUSTO,
                                          enrichr_gencode_results_GUSTO,
                                          "fig5_SC/EnrichR_Results_DGE/results_GUSTO_DGE.csv")
all_results_GUSTO.DGE <- all_results_GUSTO

## Create lollipop plots for DGE visualization
# Get list of datasets
datasets_Gen3G <- unique(all_results_Gen3G$Dataset)
datasets_GUSTO <- unique(all_results_GUSTO$Dataset)

# Create lollipop plots for each dataset
make_lolipops_perdataset(datasets_Gen3G, all_results_Gen3G,
                         "fig5_SC/EnrichR_Results_DGE/Plots/",
                         "_lollipop_DGE.pdf")

make_lolipops_perdataset(datasets_GUSTO, all_results_GUSTO,
                         "fig5_SC/EnrichR_Results_DGE/Plots/",
                         "_lollipop_DGE.pdf")

################################################################################
# Mediation Results Visualization and Functional Annotation
################################################################################

# Set up biomart for gene annotation
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

## Get top-level GO BP terms for functional categorization
bp_parents <- as.list(GOBPPARENTS)
bp_top <- as.list(GOBPCHILDREN)[["GO:0008150"]]
bp_top <- unique(bp_top)
top_bp_names <- Term(GOTERM[bp_top])
names(top_bp_names) <- bp_top

## Function to find top-level ancestor of a GO term
get_top_level_bp <- function(go_id, top_terms, parent_map) {
  visited <- character()
  while (!is.null(go_id) && !(go_id %in% visited)) {
    if (go_id %in% top_terms) return(go_id)
    visited <- c(visited, go_id)
    go_id <- parent_map[[go_id]]
    if (!is.null(go_id)) go_id <- go_id[1]  # Take first parent if multiple
  }
  return(NA)
}

# Define functional categories for visualization
category_map <- list(
  "Immune" = c("immune system process", "response to stimulus", "viral process"),
  "Development" = c("developmental process", "multicellular organismal process", "reproductive process"),
  "Regulation" = c("biological regulation", "negative regulation of biological process", 
                   "homeostatic process", "rhythmic process"),
  "Metabolism" = c("metabolic process"),
  "Cell Organization" = c("cellular process", "localization"),
  "Unknown" = c("Unknown")
)

priority <- c("Immune", "Development", "Regulation", "Metabolism", "Cell Organization", "Unknown")

# Function to assign category based on GO terms
assign_category <- function(term) {
  for (cat in names(category_map)) {
    if (term %in% category_map[[cat]]) return(cat)
  }
  return(NA)
}

################################################################################
# GUSTO Mediation Results Annotation and Visualization
################################################################################

## Process GUSTO significant mediation results
GUSTO.mediation.gdm.sig <- GUSTO.mediation.gdm[GUSTO.mediation.gdm$pvalue < 0.05, ]
GUSTO.mediation.gdm.sig <- GUSTO.mediation.gdm.sig[order(GUSTO.mediation.gdm.sig$estimate_abs, decreasing = TRUE), ]
GUSTO.mediation.gdm.sig$geneID <- sapply(strsplit(GUSTO.mediation.gdm.sig$gene_id, '[.]'), function(x) x[1])

# Map to Entrez IDs for GO analysis
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = GUSTO.mediation.gdm.sig$geneID,
                     column = "ENTREZID",
                     keytype = "ENSEMBL",
                     multiVals = "first")

# Get GO BP annotations
go_bp <- AnnotationDbi::select(org.Hs.eg.db,
                               keys = entrez_ids,
                               columns = c("GO"),
                               keytype = "ENTREZID")
go_bp <- go_bp[go_bp$ONTOLOGY == "BP", ]
go_bp$TopLevelGO <- sapply(go_bp$GO, get_top_level_bp, top_terms = bp_top, parent_map = bp_parents)
go_bp$TopLevelTerm <- NA
valid_idx <- !is.na(go_bp$TopLevelGO)
go_bp$TopLevelTerm[valid_idx] <- Term(GOTERM[go_bp$TopLevelGO[valid_idx]])
go_bp$ENSEMBL <- names(entrez_ids)[match(go_bp$ENTREZID, entrez_ids)]
go_bp <- go_bp[, c(7, 6)]
names(go_bp) <- c("geneID", "term")

go_bp[is.na(go_bp$term), "term"] <- "Unknown"
go_bp$top <- sapply(go_bp$term, assign_category)
go_bp <- go_bp[, c(1, 3)]
go_bp <- go_bp[!duplicated(go_bp), ]

# Prioritize categories per gene (take highest priority category)
go_bp <- go_bp %>%
  mutate(priority_rank = match(top, priority)) %>%
  group_by(geneID) %>%
  slice_min(order_by = priority_rank, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(geneID, top)

GUSTO.mediation.gdm.sig <- left_join(GUSTO.mediation.gdm.sig, go_bp)
names(GUSTO.mediation.gdm.sig)[11] <- "Function"
GUSTO.mediation.gdm.sig[is.na(GUSTO.mediation.gdm.sig$Function), "Function"] <- "Unknown"

# Get gene symbols for visualization
gene_symbols <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = GUSTO.mediation.gdm.sig$geneID,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

gene_symbols <- data.frame(geneID = GUSTO.mediation.gdm.sig$geneID, symbol = gene_symbols)

# Get additional gene information for genes without symbols
gene_info <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "description", "gene_biotype"),
  filters = "ensembl_gene_id",
  values = gene_symbols[is.na(gene_symbols$symbol), "geneID"],
  mart = mart
)

# Manual annotation for specific genes
gene_info[6, 2] <- "OAT-like pseudogene"
gene_info[8, 2] <- "SAFB pseudogene"
gene_info[9, 2] <- "CCNB1IP1 antisense"
gene_info[10, 2] <- "HLA-F pseudogene"
gene_info[11, 2] <- "lnc-PMM2-6"
gene_info[12, 2] <- "chr7-lncRNA"
gene_info[13, 2] <- "lnc-ZSCAN2-5"
names(gene_info)[1] <- "geneID"

gene_symbols <- left_join(gene_symbols, gene_info[, 1:2])
gene_symbols$gene_symbol <- NA
gene_symbols[is.na(gene_symbols$symbol), "gene_symbol"] <- gene_symbols[is.na(gene_symbols$symbol), "external_gene_name"]
gene_symbols[!is.na(gene_symbols$symbol), "gene_symbol"] <- gene_symbols[!is.na(gene_symbols$symbol), "symbol"]
gene_symbols <- gene_symbols[, c(1, 4)]

GUSTO.mediation.gdm.sig <- left_join(GUSTO.mediation.gdm.sig, gene_symbols)
GUSTO.mediation.gdm.sig$label <- paste0(GUSTO.mediation.gdm.sig$gene_symbol, "\n(",
                                        GUSTO.mediation.gdm.sig$transcript_id, ")")

# Prepare data for visualization
GUSTO.mediation.gdm.sig <- GUSTO.mediation.gdm.sig[order(GUSTO.mediation.gdm.sig$estimate_abs), ]
GUSTO.mediation.gdm.sig$label <- factor(GUSTO.mediation.gdm.sig$label,
                                        levels = unique(GUSTO.mediation.gdm.sig$label))

GUSTO.mediation.gdm.sig$Function <- factor(GUSTO.mediation.gdm.sig$Function,
                                           levels = c("Cell Organization", "Metabolism",
                                                      "Development", "Immune", "Regulation",
                                                      "Unknown"))

# Separate by annotation type and select top 20
GUSTO.mediation.gdm.sig.a <- GUSTO.mediation.gdm.sig[GUSTO.mediation.gdm.sig$annotation == "lr-assembly", ]
GUSTO.mediation.gdm.sig.a <- GUSTO.mediation.gdm.sig.a[1:20, ]
GUSTO.mediation.gdm.sig.a <- GUSTO.mediation.gdm.sig.a[!is.na(GUSTO.mediation.gdm.sig.a$transcript_id), ]

GUSTO.mediation.gdm.sig.g <- GUSTO.mediation.gdm.sig[GUSTO.mediation.gdm.sig$annotation == "gencode", ]
GUSTO.mediation.gdm.sig.g <- GUSTO.mediation.gdm.sig.g[1:20, ]
GUSTO.mediation.gdm.sig.g <- GUSTO.mediation.gdm.sig.g[!is.na(GUSTO.mediation.gdm.sig.g$transcript_id), ]

## Figure 4C: Top transcript mediators of GDM effects on birth weight in GUSTO (lr-assembly)
fig4c <- ggplot(GUSTO.mediation.gdm.sig.a, aes(x = estimate, y = label, color = Function)) + 
  geom_errorbarh(aes(xmax = upper_95CI, xmin = lower_95CI), size = .5, height = .2) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = c("#f7766d", "#7dad07", "#cc7dfc", "#04bcc4", "#e9a4ae", "grey")) +
  theme(panel.grid.minor = element_blank()) +
  ylab("") +
  xlab("Total indirect effect") +
  ggtitle("lr-assembly")

fig4c
ggsave("fig4c.png", fig4c, width = 5, height = 6, units = "in", dpi = 300)

################################################################################
# Gen3G Mediation Results Annotation and Visualization
################################################################################

## Process Gen3G significant mediation results
Gen3G.mediation.gdm.sig <- Gen3G.mediation.gdm[Gen3G.mediation.gdm$pvalue < 0.05, ]
Gen3G.mediation.gdm.sig <- Gen3G.mediation.gdm.sig[order(Gen3G.mediation.gdm.sig$estimate_abs, decreasing = TRUE), ]
Gen3G.mediation.gdm.sig$geneID <- sapply(strsplit(Gen3G.mediation.gdm.sig$gene_id, '[.]'), function(x) x[1])

# Map to Entrez IDs for GO analysis
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = Gen3G.mediation.gdm.sig$geneID,
                     column = "ENTREZID",
                     keytype = "ENSEMBL",
                     multiVals = "first")

# Get GO BP annotations
go_bp <- AnnotationDbi::select(org.Hs.eg.db,
                               keys = entrez_ids,
                               columns = c("GO"),
                               keytype = "ENTREZID")
go_bp <- go_bp[go_bp$ONTOLOGY == "BP", ]
go_bp$TopLevelGO <- sapply(go_bp$GO, get_top_level_bp, top_terms = bp_top, parent_map = bp_parents)
go_bp$TopLevelTerm <- NA
valid_idx <- !is.na(go_bp$TopLevelGO)
go_bp$TopLevelTerm[valid_idx] <- Term(GOTERM[go_bp$TopLevelGO[valid_idx]])
go_bp$ENSEMBL <- names(entrez_ids)[match(go_bp$ENTREZID, entrez_ids)]
go_bp <- go_bp[, c(7, 6)]
names(go_bp) <- c("geneID", "term")

go_bp[is.na(go_bp$term), "term"] <- "Unknown"
go_bp$top <- sapply(go_bp$term, assign_category)
go_bp <- go_bp[, c(1, 3)]
go_bp <- go_bp[!duplicated(go_bp), ]

# Prioritize categories per gene
go_bp <- go_bp %>%
  mutate(priority_rank = match(top, priority)) %>%
  group_by(geneID) %>%
  slice_min(order_by = priority_rank, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(geneID, top)

Gen3G.mediation.gdm.sig <- left_join(Gen3G.mediation.gdm.sig, go_bp)
names(Gen3G.mediation.gdm.sig)[11] <- "Function"
Gen3G.mediation.gdm.sig[is.na(Gen3G.mediation.gdm.sig$Function), "Function"] <- "Unknown"

# Get gene symbols for visualization
gene_symbols <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = Gen3G.mediation.gdm.sig$geneID,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

gene_symbols <- data.frame(geneID = Gen3G.mediation.gdm.sig$geneID, symbol = gene_symbols)

# Get additional gene information for genes without symbols
gene_info <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "description", "gene_biotype"),
  filters = "ensembl_gene_id",
  values = gene_symbols[is.na(gene_symbols$symbol), "geneID"],
  mart = mart
)

# Manual annotation for specific genes
gene_info[1, 2] <- "ZNF793 antisense"
gene_info[2, 2] <- "lnc-HIVEP1-1"
gene_info[3, 2] <- "lnc-ISOC1-1"
gene_info[4, 2] <- "lnc-HAND1-1"
gene_info[6, 2] <- "lnc-ACHE-11"
names(gene_info)[1] <- "geneID"

gene_symbols <- left_join(gene_symbols, gene_info[, 1:2])
gene_symbols$gene_symbol <- NA
gene_symbols[is.na(gene_symbols$symbol), "gene_symbol"] <- gene_symbols[is.na(gene_symbols$symbol), "external_gene_name"]
gene_symbols[!is.na(gene_symbols$symbol), "gene_symbol"] <- gene_symbols[!is.na(gene_symbols$symbol), "symbol"]
gene_symbols <- gene_symbols[, c(1, 4)]

Gen3G.mediation.gdm.sig <- left_join(Gen3G.mediation.gdm.sig, gene_symbols)
Gen3G.mediation.gdm.sig$label <- paste0(Gen3G.mediation.gdm.sig$gene_symbol, "\n(",
                                        Gen3G.mediation.gdm.sig$transcript_id, ")")

# Prepare data for visualization
Gen3G.mediation.gdm.sig <- Gen3G.mediation.gdm.sig[order(Gen3G.mediation.gdm.sig$estimate_abs), ]
Gen3G.mediation.gdm.sig$label <- factor(Gen3G.mediation.gdm.sig$label,
                                        levels = unique(Gen3G.mediation.gdm.sig$label))

Gen3G.mediation.gdm.sig$Function <- factor(Gen3G.mediation.gdm.sig$Function,
                                           levels = c("Cell Organization", "Metabolism",
                                                      "Development", "Immune", "Regulation",
                                                      "Unknown"))

# Separate by annotation type and select top 20
Gen3G.mediation.gdm.sig.a <- Gen3G.mediation.gdm.sig[Gen3G.mediation.gdm.sig$annotation == "lr-assembly", ]
Gen3G.mediation.gdm.sig.a <- Gen3G.mediation.gdm.sig.a[1:20, ]
Gen3G.mediation.gdm.sig.a <- Gen3G.mediation.gdm.sig.a[!is.na(Gen3G.mediation.gdm.sig.a$transcript_id), ]

Gen3G.mediation.gdm.sig.g <- Gen3G.mediation.gdm.sig[Gen3G.mediation.gdm.sig$annotation == "gencode", ]
Gen3G.mediation.gdm.sig.g <- Gen3G.mediation.gdm.sig.g[1:20, ]
Gen3G.mediation.gdm.sig.g <- Gen3G.mediation.gdm.sig.g[!is.na(Gen3G.mediation.gdm.sig.g$transcript_id), ]

## Figure 4D: Top transcript mediators of GDM effects on birth weight in Gen3G (lr-assembly)
fig4d <- ggplot(Gen3G.mediation.gdm.sig.a, aes(x = estimate, y = label, color = Function)) + 
  geom_errorbarh(aes(xmax = upper_95CI, xmin = lower_95CI), size = .5, height = .2) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = c("#f7766d", "#7dad07", "#cc7dfc", "#e9a4ae", "grey")) +
  theme(panel.grid.minor = element_blank()) +
  ylab("") +
  xlab("Total indirect effect") +
  ggtitle("lr-assembly")

fig4d
ggsave("fig4d.png", fig4d, width = 5, height = 6, units = "in", dpi = 300)

################################################################################
# CSH1 Transcript Structure Visualization (Figure 4E)
################################################################################

# Import assembly and reference annotations for transcript structure visualization
gtf.assembly <- data.frame(import("ESPRESSO_corrected_SC_filtered/ESPRESSO_corrected_SC_filtered.gtf"))
gtf.assembly$gene_id <- NULL
gtf.assembly <- left_join(gtf.assembly, tx2g.assembly)
gtf.assembly <- left_join(gtf.assembly, classif[, c("transcript_id", "structural_category")])
names(gtf.assembly)[12] <- "class"

gtf.gencode <- data.frame(import("/Users/stbresnahan/Desktop/Placenta_LRRNAseq/gencode.v45.annotation.gtf"))
gtf.gencode <- gtf.gencode[, names(gtf.gencode) %in% names(gtf.assembly)]
gtf.gencode$class <- "Reference"
gtf.gencode <- gtf.gencode[!gtf.gencode$transcript_id %in% gtf.assembly$transcript_id, ]

gtf.assembly <- rbind(gtf.assembly, gtf.gencode)
gtf.assembly$gene_id <- sapply(strsplit(gtf.assembly$gene_id, '[.]'), function(x) x[1])

# Function to create transcript structure plots
tx.plot <- function(gtf.df, gene.name, title) {
  dat <- gtf.df[gtf.df$gene_id == gene.name & gtf.df$type == "exon", ]
  dat$transcript_id <- factor(dat$transcript_id,
                              levels = unique(dat$transcript_id))
  dat$class <- factor(as.character(dat$class),
                      levels = c("Reference", "FSM", "ISM", "NIC", "NNC",
                                 "Genic Genomic", "Antisense", "Fusion",
                                 "Intergenic", "Genic Intron"))
  dat <- dat %>%
    mutate(transcript_id = factor(transcript_id,
                                  levels = dat %>%
                                    distinct(transcript_id, class) %>%
                                    arrange(class) %>%
                                    pull(transcript_id)))
  
  p <- dat %>% ggplot(aes(xstart = start, xend = end, y = transcript_id)) +
    geom_range(aes(fill = class)) +
    geom_intron(data = to_intron(dat, "transcript_id"),
                aes(strand = strand),
                arrow = grid::arrow(length = unit(0.05, "inches"))) +
    scale_y_discrete(limits = rev) +
    theme_minimal() +
    theme(panel.grid.major.y = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "top",
          axis.text.y = element_text(size = 4)) +
    ggtitle(title)
  
  return(p)
}

# Create CSH1 transcript structure plot
gtf.assembly.CSH1ir <- gtf.assembly[gtf.assembly$transcript_id %in% c("ENST00000453363.7",
                                                                      "ENST00000316193.13",
                                                                      "ENST00000329882.8",
                                                                      "ENST00000558284.1",
                                                                      "ENST00000610991.1",
                                                                      "ENST00000558661.1",
                                                                      "ESPRESSO:chr17:12985:70",
                                                                      "ESPRESSO:chr17:12985:62"), ]

## Figure 4E: CSH1 isoforms mediating GDM effects on birth weight
fig4e <- tx.plot(gtf.assembly.CSH1ir, "ENSG00000136488", 
                 "Isoforms of CSH1 (ENSG00000136488)\nmediating effect of GDM on birth weight") +
  scale_fill_manual(values = c("#66c2a5", "#e78ac3"))

fig4e
ggsave("fig4e.png", fig4e, width = 6, height = 3, units = "in", dpi = 300)

################################################################################
# End of Analysis
################################################################################