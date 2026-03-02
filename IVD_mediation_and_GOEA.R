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
# Input Datasets:
#   - Supplementary/tx2g_ESPRESSO_assembly_SC_filtered.RData
#       Transcript-to-gene mapping objects (tx2g.assembly, tx2g.gencode, tx2g.combined)
#   - Supplementary/metadata_GUSTO.RData
#       GUSTO cohort metadata (metadata.GUSTO, cordat.GUSTO)
#   - Supplementary/metadata_Gen3G.RData
#       Gen3G cohort metadata (metadata.Gen3G, cordat.Gen3G.Gen3G)
#   - Results/GUSTO_ESPRESSO_assembly_filtered_SC.RData
#       GUSTO differential expression results (all annotations and tissues)
#   - Results/Gen3G_ESPRESSO_assembly_filtered_SC.RData
#       Gen3G differential expression results (all annotations)
#   - Assembly/ESPRESSO_corrected_SC_filtered_classification.txt
#       SQANTI3 structural classification for filtered assembly
#   - Supplementary/pfam.domtblout.txt
#       HMMER protein domain predictions (Pfam database)
#   - Assembly/ESPRESSO_corrected_SC_filtered.gtf.gz
#       GTF annotation for filtered assembly (for transcript structure visualization)
#
# Output Files:
#   Data Files - Mediation Analysis:
#     - Supplementary/all_results_tables.RData
#         Comprehensive results including DE tables, mediation results, and VST data
#     - Supplementary/mediation_union_transcripts.RData
#         Mediation results for union of top mediators across cohorts
#
#   Data Files - PC Mediation Sensitivity:
#     - Supplementary/PCsensitivity.GUSTO.csv
#         PC sensitivity indirect effects (GUSTO lr-assembly)
#     - Supplementary/PCsensitivity.totals.GUSTO.csv
#         PC sensitivity total effects (GUSTO lr-assembly)
#     - Supplementary/PCsensitivity.GUSTO.gencode.csv
#         PC sensitivity indirect effects (GUSTO GENCODEv45)
#     - Supplementary/PCsensitivity.totals.GUSTO.gencode.csv
#         PC sensitivity total effects (GUSTO GENCODEv45)
#     - Supplementary/PCsensitivity.GUSTO.gencode_plus.csv
#         PC sensitivity indirect effects (GUSTO GENCODE+)
#     - Supplementary/PCsensitivity.totals.GUSTO.gencode_plus.csv
#         PC sensitivity total effects (GUSTO GENCODE+)
#     - Supplementary/PCsensitivity.GUSTO.adipose.csv
#         PC sensitivity indirect effects (GUSTO Adipose)
#     - Supplementary/PCsensitivity.totals.GUSTO.adipose.csv
#         PC sensitivity total effects (GUSTO Adipose)
#     - Supplementary/PCsensitivity.GUSTO.fibroblasts.csv
#         PC sensitivity indirect effects (GUSTO Fibroblasts)
#     - Supplementary/PCsensitivity.totals.GUSTO.fibroblasts.csv
#         PC sensitivity total effects (GUSTO Fibroblasts)
#     - Supplementary/PCsensitivity.Gen3G.csv
#         PC sensitivity indirect effects (Gen3G lr-assembly)
#     - Supplementary/PCsensitivity.totals.Gen3G.csv
#         PC sensitivity total effects (Gen3G lr-assembly)
#     - Supplementary/PCsensitivity.Gen3G.gencode.csv
#         PC sensitivity indirect effects (Gen3G GENCODEv45)
#     - Supplementary/PCsensitivity.totals.Gen3G.gencode.csv
#         PC sensitivity total effects (Gen3G GENCODEv45)
#     - Supplementary/PCsensitivity.Gen3G.gencode_plus.csv
#         PC sensitivity indirect effects (Gen3G GENCODE+)
#     - Supplementary/PCsensitivity.totals.Gen3G.gencode_plus.csv
#         PC sensitivity total effects (Gen3G GENCODE+)
#
#   Data Files - Enrichment Analysis:
#     - EnrichR_Results/EnrichR_Results_top50isorichgenes/results.csv
#         Enrichment results for top 50 isoform-rich genes
#     - top50isorichgenes_GOMF.csv
#         GO Molecular Function results for top 50 isoform-rich genes
#     - EnrichR_Results/EnrichR_Results_DTE/results_GUSTO_DTE.csv
#         GUSTO differential transcript expression enrichment results
#     - EnrichR_Results/EnrichR_Results_DTE/results_Gen3G_DTE.csv
#         Gen3G differential transcript expression enrichment results
#     - EnrichR_Results/EnrichR_Results_DGE/results_GUSTO_DGE.csv
#         GUSTO differential gene expression enrichment results
#     - EnrichR_Results/EnrichR_Results_DGE/results_Gen3G_DGE.csv
#         Gen3G differential gene expression enrichment results
#     - EnrichR_Results/EnrichR_Results_DTE/*.csv
#         Individual database results for DTE analysis
#     - EnrichR_Results/EnrichR_Results_DGE/*.csv
#         Individual database results for DEG analysis
#
#   Figures - Mediation Analysis:
#     - Figures/fig3f.pdf
#         Direct and indirect effects across cohorts and annotations (GUSTO only)
#     - Figures/fig4a.pdf
#         Mediation effects for union of top 10 significant mediators across cohorts
#     - Figures/fig4b.pdf
#         CSH1 isoform structure with protein domains (highlighting mediators)
#
#   Figures - PC Sensitivity Analysis:
#     - Figures/figs8a.pdf
#         PC mediation sensitivity for GUSTO (GENCODEv45, GENCODE+, lr-assembly)
#     - Figures/figs8b.pdf
#         PC mediation sensitivity for Gen3G (GENCODEv45, GENCODE+, lr-assembly)
#     - Figures/figs8d.pdf
#         PC mediation sensitivity for GUSTO tissues (Placenta, Adipose, Fibroblasts)
#
#   Figures - Enrichment Analysis:
#     - EnrichR_Results/EnrichR_Results_DTE/Plots/*_lollipop_DTE.pdf
#         Lollipop plots showing top enriched terms for DTE (GUSTO and Gen3G, assembly and gencode)
#     - EnrichR_Results/EnrichR_Results_DGE/Plots/*_lollipop_DGE.pdf
#         Lollipop plots showing top enriched terms for DEG (GUSTO and Gen3G, assembly and gencode)
#
################################################################################

################################################################################
# Environment Setup
################################################################################

# Configure R library paths for HPC environment
# .libPaths(c("/home/stbresnahan/R/ubuntu/4.3.1", .libPaths()))

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
load("Supplementary/tx2g_ESPRESSO_assembly_SC_filtered.RData")

# Load short-read cohort metadata
load("Supplementary/metadata_GUSTO.RData")
load("Supplementary/metadata_Gen3G.RData")

# Load differential expression results from both cohorts
load("Results/GUSTO_ESPRESSO_assembly_filtered_SC.RData")
load("Results/Gen3G_ESPRESSO_assembly_filtered_SC.RData")

# Import assembly isoform classification from SQANTI3
classif <- read.table("Assembly/ESPRESSO_corrected_SC_filtered_classification.txt", 
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
GUSTO.vst.gencode_plus <- assay(varianceStabilizingTransformation(GUSTO.drgC.gencode_plus.gdm[[3]], blind = FALSE))
Gen3G.vst.assembly <- assay(varianceStabilizingTransformation(Gen3G.drgC.assembly.gdm[[3]], blind = FALSE))
Gen3G.vst.gencode <- assay(varianceStabilizingTransformation(Gen3G.drgC.gencode.gdm[[3]], blind = FALSE))
Gen3G.vst.gencode_plus <- assay(varianceStabilizingTransformation(Gen3G.drgC.gencode_plus.gdm[[3]], blind = FALSE))

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

results.gencode_plus.GUSTO.gdm <- GUSTO.drgC.gencode_plus.gdm[[4]]
results.gencode_plus.GUSTO.gdm$transcript_id <- row.names(results.gencode_plus.GUSTO.gdm)
results.gencode_plus.Gen3G.gdm <- Gen3G.drgC.gencode_plus.gdm[[4]]
results.gencode_plus.Gen3G.gdm$transcript_id <- row.names(results.gencode_plus.Gen3G.gdm)

# Extract gene-level aggregated results for transcript-level analysis
results.assembly.GUSTO.gdm.gene <- GUSTO.drgC.assembly.gdm[[5]]
results.assembly.Gen3G.gdm.gene <- Gen3G.drgC.assembly.gdm[[5]]
results.gencode.GUSTO.gdm.gene <- GUSTO.drgC.gencode.gdm[[5]]
results.gencode.Gen3G.gdm.gene <- Gen3G.drgC.gencode.gdm[[5]]
results.gencode_plus.GUSTO.gdm.gene <- GUSTO.drgC.gencode_plus.gdm[[5]]
results.gencode_plus.Gen3G.gdm.gene <- Gen3G.drgC.gencode_plus.gdm[[5]]

# Extract gene-level differential expression results
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

# Define significant gene/transcript lists at FDR < 0.1
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
# Individual Transcript Mediation Analysis
################################################################################

# Set up parallel backend for computationally intensive mediation analysis
num_cores <- detectCores() - 1  # Leave one core free for system processes
cl <- makeCluster(num_cores)
registerDoParallel(cl)

################################################################################
# GUSTO
################################################################################

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
GUSTO.expr.assembly$birth_weight <- as.numeric(GUSTO.expr.assembly$birth_weight)
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
                                                                  boot = TRUE, sims = 100)
                                        
                                        # Extract coefficients
                                        med_summary <- summary(mediator_model)
                                        out_summary <- summary(outcome_model)
                                        
                                        a_est <- coef(med_summary)["gdm1", "Estimate"]
                                        a_se <- coef(med_summary)["gdm1", "Std. Error"]
                                        a_pval <- coef(med_summary)["gdm1", "Pr(>|t|)"]
                                        
                                        b_est <- coef(out_summary)["TXE", "Estimate"]
                                        b_se <- coef(out_summary)["TXE", "Std. Error"]
                                        b_pval <- coef(out_summary)["TXE", "Pr(>|t|)"]
                                        
                                        data.frame(transcript_id = txids[tx],
                                                   estimate = mediate_result$d0,
                                                   lower_95CI = mediate_result$d0.ci[1],
                                                   upper_95CI = mediate_result$d0.ci[2],
                                                   pvalue = mediate_result$d0.p,
                                                   a_path_estimate = a_est,
                                                   a_path_lower_95CI = a_est - 1.96 * a_se,
                                                   a_path_upper_95CI = a_est + 1.96 * a_se,
                                                   a_path_pvalue = a_pval,
                                                   b_path_estimate = b_est,
                                                   b_path_lower_95CI = b_est - 1.96 * b_se,
                                                   b_path_upper_95CI = b_est + 1.96 * b_se,
                                                   b_path_pvalue = b_pval)
                                      }, error = function(e) {
                                        warning(paste("Error in transcript", txids[tx], ":", e$message))
                                        data.frame(transcript_id = txids[tx],
                                                   estimate = NA,
                                                   lower_95CI = NA,
                                                   upper_95CI = NA,
                                                   pvalue = NA,
                                                   a_path_estimate = NA,
                                                   a_path_lower_95CI = NA,
                                                   a_path_upper_95CI = NA,
                                                   a_path_pvalue = NA,
                                                   b_path_estimate = NA,
                                                   b_path_lower_95CI = NA,
                                                   b_path_upper_95CI = NA,
                                                   b_path_pvalue = NA)
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
GUSTO.expr.gencode$birth_weight <- as.numeric(GUSTO.expr.gencode$birth_weight)
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
                                                                 boot = TRUE, sims = 100)
                                       
                                       # Extract coefficients
                                       med_summary <- summary(mediator_model)
                                       out_summary <- summary(outcome_model)
                                       
                                       a_est <- coef(med_summary)["gdm1", "Estimate"]
                                       a_se <- coef(med_summary)["gdm1", "Std. Error"]
                                       a_pval <- coef(med_summary)["gdm1", "Pr(>|t|)"]
                                       
                                       b_est <- coef(out_summary)["TXE", "Estimate"]
                                       b_se <- coef(out_summary)["TXE", "Std. Error"]
                                       b_pval <- coef(out_summary)["TXE", "Pr(>|t|)"]
                                       
                                       data.frame(transcript_id = txids[tx],
                                                  estimate = mediate_result$d0,
                                                  lower_95CI = mediate_result$d0.ci[1],
                                                  upper_95CI = mediate_result$d0.ci[2],
                                                  pvalue = mediate_result$d0.p,
                                                  a_path_estimate = a_est,
                                                  a_path_lower_95CI = a_est - 1.96 * a_se,
                                                  a_path_upper_95CI = a_est + 1.96 * a_se,
                                                  a_path_pvalue = a_pval,
                                                  b_path_estimate = b_est,
                                                  b_path_lower_95CI = b_est - 1.96 * b_se,
                                                  b_path_upper_95CI = b_est + 1.96 * b_se,
                                                  b_path_pvalue = b_pval)
                                     }, error = function(e) {
                                       warning(paste("Error in transcript", txids[tx], ":", e$message))
                                       data.frame(transcript_id = txids[tx],
                                                  estimate = NA,
                                                  lower_95CI = NA,
                                                  upper_95CI = NA,
                                                  pvalue = NA,
                                                  a_path_estimate = NA,
                                                  a_path_lower_95CI = NA,
                                                  a_path_upper_95CI = NA,
                                                  a_path_pvalue = NA,
                                                  b_path_estimate = NA,
                                                  b_path_lower_95CI = NA,
                                                  b_path_upper_95CI = NA,
                                                  b_path_pvalue = NA)
                                     })
                                   }

row.names(GUSTO.mediation.gencode) <- NULL
GUSTO.mediation.gencode$estimate_abs <- abs(GUSTO.mediation.gencode$estimate)
GUSTO.mediation.gencode <- left_join(GUSTO.mediation.gencode, tx2g.gencode)


## GUSTO gencode_plus: Individual transcript mediation analysis
GUSTO.expr.gencode_plus <- data.frame(GUSTO.vst.gencode_plus)
GUSTO.expr.gencode_plus$transcript_id <- row.names(GUSTO.expr.gencode_plus)

# Filter to significant differentially expressed transcripts
GUSTO.expr.gencode_plus <- GUSTO.expr.gencode_plus[GUSTO.expr.gencode_plus$transcript_id %in% 
                                                     results.gencode_plus.GUSTO.gdm[results.gencode_plus.GUSTO.gdm$padj < 0.1, "transcript_id"], ]

# Reshape data for mediation analysis
GUSTO.expr.gencode_plus <- GUSTO.expr.gencode_plus %>%
  pivot_longer(
    cols = -transcript_id,
    names_to = "SampleID",
    values_to = "TXE")

# Add covariate data
GUSTO.expr.gencode_plus <- left_join(GUSTO.expr.gencode_plus, 
                                     cordat.GUSTO[, c("SampleID", "gdm", "birth_weight", "sex", "GA")])
GUSTO.expr.gencode_plus$birth_weight <- as.numeric(GUSTO.expr.gencode_plus$birth_weight)
txids <- unique(GUSTO.expr.gencode_plus$transcript_id)

# Parallel mediation analysis for each transcript
GUSTO.mediation.gencode_plus <- foreach(tx = 1:length(txids), 
                                        .packages = c("mediation"),
                                        .combine = rbind) %dopar% {
                                          
                                          data <- GUSTO.expr.gencode_plus[GUSTO.expr.gencode_plus$transcript_id == txids[tx], ]
                                          mediator_model <- lm(TXE ~ gdm + sex + GA, data = data)
                                          outcome_model <- lm(birth_weight ~ gdm + TXE + sex + GA, data = data)
                                          
                                          tryCatch({
                                            mediate_result <- mediate(mediator_model, outcome_model, 
                                                                      treat = "gdm", mediator = "TXE", 
                                                                      boot = TRUE, sims = 100)
                                            
                                            # Extract coefficients
                                            med_summary <- summary(mediator_model)
                                            out_summary <- summary(outcome_model)
                                            
                                            a_est <- coef(med_summary)["gdm1", "Estimate"]
                                            a_se <- coef(med_summary)["gdm1", "Std. Error"]
                                            a_pval <- coef(med_summary)["gdm1", "Pr(>|t|)"]
                                            
                                            b_est <- coef(out_summary)["TXE", "Estimate"]
                                            b_se <- coef(out_summary)["TXE", "Std. Error"]
                                            b_pval <- coef(out_summary)["TXE", "Pr(>|t|)"]
                                            
                                            data.frame(transcript_id = txids[tx],
                                                       estimate = mediate_result$d0,
                                                       lower_95CI = mediate_result$d0.ci[1],
                                                       upper_95CI = mediate_result$d0.ci[2],
                                                       pvalue = mediate_result$d0.p,
                                                       a_path_estimate = a_est,
                                                       a_path_lower_95CI = a_est - 1.96 * a_se,
                                                       a_path_upper_95CI = a_est + 1.96 * a_se,
                                                       a_path_pvalue = a_pval,
                                                       b_path_estimate = b_est,
                                                       b_path_lower_95CI = b_est - 1.96 * b_se,
                                                       b_path_upper_95CI = b_est + 1.96 * b_se,
                                                       b_path_pvalue = b_pval)
                                          }, error = function(e) {
                                            warning(paste("Error in transcript", txids[tx], ":", e$message))
                                            data.frame(transcript_id = txids[tx],
                                                       estimate = NA,
                                                       lower_95CI = NA,
                                                       upper_95CI = NA,
                                                       pvalue = NA,
                                                       a_path_estimate = NA,
                                                       a_path_lower_95CI = NA,
                                                       a_path_upper_95CI = NA,
                                                       a_path_pvalue = NA,
                                                       b_path_estimate = NA,
                                                       b_path_lower_95CI = NA,
                                                       b_path_upper_95CI = NA,
                                                       b_path_pvalue = NA)
                                          })
                                        }

row.names(GUSTO.mediation.gencode_plus) <- NULL
GUSTO.mediation.gencode_plus$estimate_abs <- abs(GUSTO.mediation.gencode_plus$estimate)
GUSTO.mediation.gencode_plus <- left_join(GUSTO.mediation.gencode_plus, tx2g.gencode)

# Combine GUSTO mediation results
GUSTO.mediation.assembly$annotation <- "lr-assembly"
GUSTO.mediation.gencode$annotation <- "gencode"
GUSTO.mediation.gencode_plus$annotation <- "gencode_plus"

GUSTO.mediation.gdm <- rbind(GUSTO.mediation.gencode, GUSTO.mediation.gencode_plus, GUSTO.mediation.assembly)
GUSTO.mediation.gdm[GUSTO.mediation.gdm$pvalue == 0, "pvalue"] <- 2e-16
GUSTO.mediation.gdm$cohort <- "GUSTO"


################################################################################
# Additional VST for Adipose and Fibroblasts
################################################################################

GUSTO.vst.adipose <- assay(varianceStabilizingTransformation(GUSTO.drgC.adipose.gdm[[3]], blind = FALSE))
GUSTO.vst.fibroblasts <- assay(varianceStabilizingTransformation(GUSTO.drgC.fibroblasts.gdm[[3]], blind = FALSE))

################################################################################
# Additional Results Tables for Adipose and Fibroblasts
################################################################################

results.adipose.GUSTO.gdm <- GUSTO.drgC.adipose.gdm[[4]]
results.adipose.GUSTO.gdm$transcript_id <- row.names(results.adipose.GUSTO.gdm)
results.fibroblasts.GUSTO.gdm <- GUSTO.drgC.fibroblasts.gdm[[4]]
results.fibroblasts.GUSTO.gdm$transcript_id <- row.names(results.fibroblasts.GUSTO.gdm)

results.adipose.GUSTO.gdm.gene <- GUSTO.drgC.adipose.gdm[[5]]
results.fibroblasts.GUSTO.gdm.gene <- GUSTO.drgC.fibroblasts.gdm[[5]]

results.adipose.GUSTO.gdm.gene <- GUSTO.drgC.adipose.gdm.gene[[4]]
results.adipose.GUSTO.gdm.gene$gene_id <- row.names(results.adipose.GUSTO.gdm.gene)
results.fibroblasts.GUSTO.gdm.gene <- GUSTO.drgC.fibroblasts.gdm.gene[[4]]
results.fibroblasts.GUSTO.gdm.gene$gene_id <- row.names(results.fibroblasts.GUSTO.gdm.gene)

GUSTO.gdm.DTE.adipose <- results.adipose.GUSTO.gdm[results.adipose.GUSTO.gdm$padj < 0.1, "transcript_id"]
GUSTO.gdm.DTE.fibroblasts <- results.fibroblasts.GUSTO.gdm[results.fibroblasts.GUSTO.gdm$padj < 0.1, "transcript_id"]
GUSTO.gdm.DGE.adipose <- results.adipose.GUSTO.gdm.gene[results.adipose.GUSTO.gdm.gene$padj < 0.1, "gene_id"]
GUSTO.gdm.DGE.fibroblasts <- results.fibroblasts.GUSTO.gdm.gene[results.fibroblasts.GUSTO.gdm.gene$padj < 0.1, "gene_id"]

################################################################################
# Gen3G
################################################################################

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
Gen3G.expr.assembly$birth_weight <- as.numeric(Gen3G.expr.assembly$birth_weight)
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
                                        
                                        # Extract coefficients
                                        med_summary <- summary(mediator_model)
                                        out_summary <- summary(outcome_model)
                                        
                                        a_est <- coef(med_summary)["gdm1", "Estimate"]
                                        a_se <- coef(med_summary)["gdm1", "Std. Error"]
                                        a_pval <- coef(med_summary)["gdm1", "Pr(>|t|)"]
                                        
                                        b_est <- coef(out_summary)["TXE", "Estimate"]
                                        b_se <- coef(out_summary)["TXE", "Std. Error"]
                                        b_pval <- coef(out_summary)["TXE", "Pr(>|t|)"]
                                        
                                        data.frame(transcript_id = txids[tx],
                                                   estimate = mediate_result$d0,
                                                   lower_95CI = mediate_result$d0.ci[1],
                                                   upper_95CI = mediate_result$d0.ci[2],
                                                   pvalue = mediate_result$d0.p,
                                                   a_path_estimate = a_est,
                                                   a_path_lower_95CI = a_est - 1.96 * a_se,
                                                   a_path_upper_95CI = a_est + 1.96 * a_se,
                                                   a_path_pvalue = a_pval,
                                                   b_path_estimate = b_est,
                                                   b_path_lower_95CI = b_est - 1.96 * b_se,
                                                   b_path_upper_95CI = b_est + 1.96 * b_se,
                                                   b_path_pvalue = b_pval)
                                      }, error = function(e) {
                                        warning(paste("Error in transcript", txids[tx], ":", e$message))
                                        data.frame(transcript_id = txids[tx],
                                                   estimate = NA,
                                                   lower_95CI = NA,
                                                   upper_95CI = NA,
                                                   pvalue = NA,
                                                   a_path_estimate = NA,
                                                   a_path_lower_95CI = NA,
                                                   a_path_upper_95CI = NA,
                                                   a_path_pvalue = NA,
                                                   b_path_estimate = NA,
                                                   b_path_lower_95CI = NA,
                                                   b_path_upper_95CI = NA,
                                                   b_path_pvalue = NA)
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
Gen3G.expr.gencode$birth_weight <- as.numeric(Gen3G.expr.gencode$birth_weight)
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
                                       
                                       # Extract coefficients
                                       med_summary <- summary(mediator_model)
                                       out_summary <- summary(outcome_model)
                                       
                                       a_est <- coef(med_summary)["gdm1", "Estimate"]
                                       a_se <- coef(med_summary)["gdm1", "Std. Error"]
                                       a_pval <- coef(med_summary)["gdm1", "Pr(>|t|)"]
                                       
                                       b_est <- coef(out_summary)["TXE", "Estimate"]
                                       b_se <- coef(out_summary)["TXE", "Std. Error"]
                                       b_pval <- coef(out_summary)["TXE", "Pr(>|t|)"]
                                       
                                       data.frame(transcript_id = txids[tx],
                                                  estimate = mediate_result$d0,
                                                  lower_95CI = mediate_result$d0.ci[1],
                                                  upper_95CI = mediate_result$d0.ci[2],
                                                  pvalue = mediate_result$d0.p,
                                                  a_path_estimate = a_est,
                                                  a_path_lower_95CI = a_est - 1.96 * a_se,
                                                  a_path_upper_95CI = a_est + 1.96 * a_se,
                                                  a_path_pvalue = a_pval,
                                                  b_path_estimate = b_est,
                                                  b_path_lower_95CI = b_est - 1.96 * b_se,
                                                  b_path_upper_95CI = b_est + 1.96 * b_se,
                                                  b_path_pvalue = b_pval)
                                     }, error = function(e) {
                                       warning(paste("Error in transcript", txids[tx], ":", e$message))
                                       data.frame(transcript_id = txids[tx],
                                                  estimate = NA,
                                                  lower_95CI = NA,
                                                  upper_95CI = NA,
                                                  pvalue = NA,
                                                  a_path_estimate = NA,
                                                  a_path_lower_95CI = NA,
                                                  a_path_upper_95CI = NA,
                                                  a_path_pvalue = NA,
                                                  b_path_estimate = NA,
                                                  b_path_lower_95CI = NA,
                                                  b_path_upper_95CI = NA,
                                                  b_path_pvalue = NA)
                                     })
                                   }

row.names(Gen3G.mediation.gencode) <- NULL
Gen3G.mediation.gencode$estimate_abs <- abs(Gen3G.mediation.gencode$estimate)
Gen3G.mediation.gencode <- left_join(Gen3G.mediation.gencode, tx2g.gencode, by = "transcript_id")


## Gen3G gencode_plus: Individual transcript mediation analysis
Gen3G.expr.gencode_plus <- data.frame(Gen3G.vst.gencode_plus)
Gen3G.expr.gencode_plus$transcript_id <- row.names(Gen3G.expr.gencode_plus)

# Filter to significant differentially expressed transcripts
Gen3G.expr.gencode_plus <- Gen3G.expr.gencode_plus[Gen3G.expr.gencode_plus$transcript_id %in% 
                                                     results.gencode_plus.Gen3G.gdm[results.gencode_plus.Gen3G.gdm$padj < 0.1, "transcript_id"], ]

# Reshape data for mediation analysis
Gen3G.expr.gencode_plus <- Gen3G.expr.gencode_plus %>%
  pivot_longer(
    cols = -transcript_id,
    names_to = "Run",
    values_to = "TXE")

# Add covariate data
Gen3G.expr.gencode_plus <- left_join(Gen3G.expr.gencode_plus, 
                                     cordat.Gen3G[, c("Run", "gdm", "birth_weight", "GA")])
Gen3G.expr.gencode_plus$birth_weight <- as.numeric(Gen3G.expr.gencode_plus$birth_weight)
txids <- unique(Gen3G.expr.gencode_plus$transcript_id)

# Parallel mediation analysis for each transcript
Gen3G.mediation.gencode_plus <- foreach(tx = 1:length(txids), 
                                        .packages = c("mediation"),
                                        .combine = rbind) %dopar% {
                                          
                                          data <- Gen3G.expr.gencode_plus[Gen3G.expr.gencode_plus$transcript_id == txids[tx], ]
                                          mediator_model <- lm(TXE ~ gdm, data = data)
                                          outcome_model <- lm(birth_weight ~ gdm + TXE, data = data)
                                          
                                          tryCatch({
                                            mediate_result <- mediate(mediator_model, outcome_model, 
                                                                      treat = "gdm", mediator = "TXE", 
                                                                      boot = TRUE, sims = 100)
                                            
                                            # Extract coefficients
                                            med_summary <- summary(mediator_model)
                                            out_summary <- summary(outcome_model)
                                            
                                            a_est <- coef(med_summary)["gdm1", "Estimate"]
                                            a_se <- coef(med_summary)["gdm1", "Std. Error"]
                                            a_pval <- coef(med_summary)["gdm1", "Pr(>|t|)"]
                                            
                                            b_est <- coef(out_summary)["TXE", "Estimate"]
                                            b_se <- coef(out_summary)["TXE", "Std. Error"]
                                            b_pval <- coef(out_summary)["TXE", "Pr(>|t|)"]
                                            
                                            data.frame(transcript_id = txids[tx],
                                                       estimate = mediate_result$d0,
                                                       lower_95CI = mediate_result$d0.ci[1],
                                                       upper_95CI = mediate_result$d0.ci[2],
                                                       pvalue = mediate_result$d0.p,
                                                       a_path_estimate = a_est,
                                                       a_path_lower_95CI = a_est - 1.96 * a_se,
                                                       a_path_upper_95CI = a_est + 1.96 * a_se,
                                                       a_path_pvalue = a_pval,
                                                       b_path_estimate = b_est,
                                                       b_path_lower_95CI = b_est - 1.96 * b_se,
                                                       b_path_upper_95CI = b_est + 1.96 * b_se,
                                                       b_path_pvalue = b_pval)
                                          }, error = function(e) {
                                            warning(paste("Error in transcript", txids[tx], ":", e$message))
                                            data.frame(transcript_id = txids[tx],
                                                       estimate = NA,
                                                       lower_95CI = NA,
                                                       upper_95CI = NA,
                                                       pvalue = NA,
                                                       a_path_estimate = NA,
                                                       a_path_lower_95CI = NA,
                                                       a_path_upper_95CI = NA,
                                                       a_path_pvalue = NA,
                                                       b_path_estimate = NA,
                                                       b_path_lower_95CI = NA,
                                                       b_path_upper_95CI = NA,
                                                       b_path_pvalue = NA)
                                          })
                                        }

row.names(Gen3G.mediation.gencode_plus) <- NULL
Gen3G.mediation.gencode_plus$estimate_abs <- abs(Gen3G.mediation.gencode_plus$estimate)
Gen3G.mediation.gencode_plus <- left_join(Gen3G.mediation.gencode_plus, tx2g.gencode, by = "transcript_id")

# Combine Gen3G mediation results
Gen3G.mediation.assembly$annotation <- "lr-assembly"
Gen3G.mediation.gencode$annotation <- "gencode"
Gen3G.mediation.gencode_plus$annotation <- "gencode_plus"

Gen3G.mediation.gdm <- rbind(Gen3G.mediation.gencode, Gen3G.mediation.gencode_plus, Gen3G.mediation.assembly)
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
  
  results_indirect <- list()
  results_total <- list()
  
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
      dplyr::select(est = est, ci.lower = ci.lower, ci.upper = ci.upper, pvalue = pvalue)
    
    # Extract total effect and proportion mediated statistics
    total_row <- param_est %>% filter(label == "total") %>%
      dplyr::select(est = est, ci.lower = ci.lower, ci.upper = ci.upper, pvalue = pvalue)
    
    prop_row <- param_est %>% filter(label == "prop_mediated") %>%
      dplyr::select(est = est, ci.lower = ci.lower, ci.upper = ci.upper, pvalue = pvalue)
    
    # Save indirect effect results
    results_indirect[[i]] <- data.frame(
      PCs = i,
      percent_variance = round(pca_cumvar[i] * 100, 2),
      total_indirect = ind_row$est,
      CI_lower = ind_row$ci.lower,
      CI_upper = ind_row$ci.upper,
      p = ind_row$pvalue
    )
    
    # Save total effect results
    results_total[[i]] <- data.frame(
      PCs = i,
      percent_variance = round(pca_cumvar[i] * 100, 2),
      total_effect = total_row$est,
      CI_lower = total_row$ci.lower,
      CI_upper = total_row$ci.upper,
      p = total_row$pvalue,
      prop_mediated = prop_row$est,
      prop_CI_lower = prop_row$ci.lower,
      prop_CI_upper = prop_row$ci.upper,
      prop_p = prop_row$pvalue
    )
  }
  
  final_results <- list(
    indirect = bind_rows(results_indirect),
    total = bind_rows(results_total)
  )
  
  return(final_results)
}

################################################################################
# GUSTO Principal Component Mediation Analysis
################################################################################

## GUSTO assembly: PC-based mediation analysis
GUSTO.expr.assembly <- data.frame(GUSTO.vst.assembly)
GUSTO.expr.assembly$transcript_id <- row.names(GUSTO.expr.assembly)

# Filter to significant transcripts for PC analysis
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
GUSTO.assembly.PCsensitivity.all <- PCmed.sensitivity(GUSTO.expr.assembly, cordat.GUSTO, 10)
GUSTO.assembly.PCsensitivity.totals <- GUSTO.assembly.PCsensitivity.all$total
GUSTO.assembly.PCsensitivity <- GUSTO.assembly.PCsensitivity.all$indirect

write.csv(GUSTO.assembly.PCsensitivity.totals,
          "Supplementary/PCsensitivity.totals.GUSTO.csv",row.names=F)
write.csv(GUSTO.assembly.PCsensitivity,
          "Supplementary/PCsensitivity.GUSTO.csv",row.names=F)


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
GUSTO.gencode.PCsensitivity.all <- PCmed.sensitivity(GUSTO.expr.gencode, cordat.GUSTO, 10)
GUSTO.gencode.PCsensitivity.totals <- GUSTO.gencode.PCsensitivity.all$total
GUSTO.gencode.PCsensitivity <- GUSTO.gencode.PCsensitivity.all$indirect

write.csv(GUSTO.gencode.PCsensitivity.totals,
          "Supplementary/PCsensitivity.totals.GUSTO.gencode.csv",row.names=F)
write.csv(GUSTO.gencode.PCsensitivity,
          "Supplementary/PCsensitivity.GUSTO.gencode.csv",row.names=F)

## GUSTO gencode_plus: PC-based mediation analysis
GUSTO.expr.gencode_plus <- data.frame(GUSTO.vst.gencode_plus)
GUSTO.expr.gencode_plus$transcript_id <- row.names(GUSTO.expr.gencode_plus)

# Filter to most significant transcripts for PC analysis
GUSTO.expr.gencode_plus <- GUSTO.expr.gencode_plus[GUSTO.expr.gencode_plus$transcript_id %in% 
                                                     results.gencode_plus.GUSTO.gdm[results.gencode_plus.GUSTO.gdm$pvalue < 0.00005, "transcript_id"], ]
GUSTO.expr.gencode_plus$transcript_id <- NULL

## Examine principal components 1-10
pca_result <- prcomp(t(as.matrix(GUSTO.expr.gencode_plus)))
variance_explained <- pca_result$sdev^2
percent_variance <- variance_explained / sum(variance_explained) * 100
percent_variance[1:10]
cumulative_variance <- cumsum(percent_variance)
cumulative_variance[1:10]
cat("PC1-7 individual variance explained:\n")
round(percent_variance[1:7], 2)
cat("\nPC1-7 cumulative variance explained:\n")
round(cumulative_variance[1:7], 2)
cat("\nTotal variance explained by PC1-7:", round(cumulative_variance[7], 2), "%\n")

plot(percent_variance[1:20], type = "b", 
     xlab = "Principal Component", 
     ylab = "Percentage Variance Explained")

## Examine correlation with maternal ethnicity
pca_result <- prcomp(t(as.matrix(GUSTO.expr.gencode_plus)))
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

# Run PC sensitivity analysis for GUSTO gencode_plus
GUSTO.gencode_plus.PCsensitivity.all <- PCmed.sensitivity(GUSTO.expr.gencode_plus, cordat.GUSTO, 10)
GUSTO.gencode_plus.PCsensitivity.totals <- GUSTO.gencode_plus.PCsensitivity.all$total
GUSTO.gencode_plus.PCsensitivity <- GUSTO.gencode_plus.PCsensitivity.all$indirect

write.csv(GUSTO.gencode_plus.PCsensitivity.totals,
          "Supplementary/PCsensitivity.totals.GUSTO.gencode_plus.csv",row.names=F)
write.csv(GUSTO.gencode_plus.PCsensitivity,
          "Supplementary/PCsensitivity.GUSTO.gencode_plus.csv",row.names=F)


## GUSTO adipose: PC-based mediation analysis
GUSTO.expr.adipose <- data.frame(GUSTO.vst.adipose)
GUSTO.expr.adipose$transcript_id <- row.names(GUSTO.expr.adipose)

# Filter to most significant transcripts for PC analysis
GUSTO.expr.adipose <- GUSTO.expr.adipose[GUSTO.expr.adipose$transcript_id %in% 
                                           results.adipose.GUSTO.gdm[results.adipose.GUSTO.gdm$pvalue < 0.00005, "transcript_id"], ]
GUSTO.expr.adipose$transcript_id <- NULL

## Examine principal components 1-10
pca_result <- prcomp(t(as.matrix(GUSTO.expr.adipose)))
variance_explained <- pca_result$sdev^2
percent_variance <- variance_explained / sum(variance_explained) * 100
percent_variance[1:10]
cumulative_variance <- cumsum(percent_variance)
cumulative_variance[1:10]
cat("PC1-10 individual variance explained:\n")
round(percent_variance[1:10], 2)
cat("\nPC1-10 cumulative variance explained:\n")
round(cumulative_variance[1:10], 2)
cat("\nTotal variance explained by PC1-10:", round(cumulative_variance[10], 2), "%\n")

plot(percent_variance[1:10], type = "b", 
     xlab = "Principal Component", 
     ylab = "Percentage Variance Explained")

## Examine correlation with maternal ethnicity (potential confounder)
pca_result <- prcomp(t(as.matrix(GUSTO.expr.adipose)))
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

# Run PC sensitivity analysis for GUSTO adipose (use 1 PC for comparison)
GUSTO.adipose.PCsensitivity.all <- PCmed.sensitivity(GUSTO.expr.adipose, cordat.GUSTO, 10)
GUSTO.adipose.PCsensitivity.totals <- GUSTO.adipose.PCsensitivity.all$total
GUSTO.adipose.PCsensitivity <- GUSTO.adipose.PCsensitivity.all$indirect

write.csv(GUSTO.adipose.PCsensitivity.totals,
          "Supplementary/PCsensitivity.totals.GUSTO.adipose.csv",row.names=F)
write.csv(GUSTO.adipose.PCsensitivity,
          "Supplementary/PCsensitivity.GUSTO.adipose.csv",row.names=F)

## GUSTO fibroblasts: PC-based mediation analysis
GUSTO.expr.fibroblasts <- data.frame(GUSTO.vst.fibroblasts)
GUSTO.expr.fibroblasts$transcript_id <- row.names(GUSTO.expr.fibroblasts)

# Filter to most significant transcripts for PC analysis
GUSTO.expr.fibroblasts <- GUSTO.expr.fibroblasts[GUSTO.expr.fibroblasts$transcript_id %in% 
                                                   results.fibroblasts.GUSTO.gdm[results.fibroblasts.GUSTO.gdm$pvalue < 0.00005, "transcript_id"], ]
GUSTO.expr.fibroblasts$transcript_id <- NULL

## Examine principal components 1-10
pca_result <- prcomp(t(as.matrix(GUSTO.expr.fibroblasts)))
variance_explained <- pca_result$sdev^2
percent_variance <- variance_explained / sum(variance_explained) * 100
percent_variance[1:10]
cumulative_variance <- cumsum(percent_variance)
cumulative_variance[1:10]
cat("PC1-10 individual variance explained:\n")
round(percent_variance[1:10], 2)
cat("\nPC1-10 cumulative variance explained:\n")
round(cumulative_variance[1:10], 2)
cat("\nTotal variance explained by PC1-10:", round(cumulative_variance[10], 2), "%\n")

plot(percent_variance[1:10], type = "b", 
     xlab = "Principal Component", 
     ylab = "Percentage Variance Explained")

## Examine correlation with maternal ethnicity
pca_result <- prcomp(t(as.matrix(GUSTO.expr.fibroblasts)))
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

# Run PC sensitivity analysis for GUSTO fibroblasts (use 2 PCs for comparison)
GUSTO.fibroblasts.PCsensitivity.all <- PCmed.sensitivity(GUSTO.expr.fibroblasts, cordat.GUSTO, 10)
GUSTO.fibroblasts.PCsensitivity.totals <- GUSTO.fibroblasts.PCsensitivity.all$total
GUSTO.fibroblasts.PCsensitivity <- GUSTO.fibroblasts.PCsensitivity.all$indirect

write.csv(GUSTO.fibroblasts.PCsensitivity.totals,
          "Supplementary/PCsensitivity.totals.GUSTO.fibroblasts.csv",row.names=F)
write.csv(GUSTO.fibroblasts.PCsensitivity,
          "Supplementary/PCsensitivity.GUSTO.fibroblasts.csv",row.names=F)




## Create sensitivity plot for GUSTO including all tissues
GUSTO.assembly.PCsensitivity$Annotation <- "Placenta (lr-assembly)"
GUSTO.gencode.PCsensitivity$Annotation <- "GENCODEv45"
GUSTO.gencode_plus.PCsensitivity$Annotation <- "GENCODE+"
GUSTO.adipose.PCsensitivity$Annotation <- "Adipose - Subcutaneous"
GUSTO.fibroblasts.PCsensitivity$Annotation <- "Cells - Cultured Fibroblasts"

PCsensitivity <- rbind(GUSTO.assembly.PCsensitivity,
                       GUSTO.gencode.PCsensitivity,
                       GUSTO.gencode_plus.PCsensitivity,
                       GUSTO.adipose.PCsensitivity,
                       GUSTO.fibroblasts.PCsensitivity)
rownames(PCsensitivity) <- NULL

PCsensitivity$Annotation <- factor(PCsensitivity$Annotation,
                                   levels = c("GENCODEv45", "GENCODE+", "Placenta (lr-assembly)",
                                              "Cells - Cultured Fibroblasts", "Adipose - Subcutaneous"))


################################################################################
# Figure S8A - PC mediation sensitivity analysis for GUSTO (main analysis)
################################################################################

PCsensitivity_a <- PCsensitivity %>%
  filter(Annotation %in% c("GENCODEv45", "GENCODE+", "Placenta (lr-assembly)"))

# Prepare data with estimate as category
PCsensitivity_a_long <- PCsensitivity_a %>%
  pivot_longer(cols = c(total_indirect, percent_variance),
               names_to = "estimate",
               values_to = "value") %>%
  mutate(estimate = factor(estimate, 
                           levels = c("total_indirect", "percent_variance"),
                           labels = c("Indirect Effect", "% Variance Explained")))

# Adjust CI columns for variance (set to NA since variance doesn't have CIs)
PCsensitivity_a_long <- PCsensitivity_a_long %>%
  mutate(CI_lower = ifelse(estimate == "Indirect Effect", CI_lower, NA_real_),
         CI_upper = ifelse(estimate == "Indirect Effect", CI_upper, NA_real_))

# Create data for conditional horizontal lines
hline_data <- data.frame(
  estimate = factor(c("Indirect Effect", "% Variance Explained"),
                    levels = c("Indirect Effect", "% Variance Explained")),
  yintercept = c(0, 75)
)

figs8a <- ggplot(PCsensitivity_a_long, aes(x = PCs, y = value, color = Annotation, fill = Annotation)) +
  geom_hline(data = hline_data, aes(yintercept = yintercept), linetype = "dotted") + 
  geom_point(size = 1.5, position = position_dodge(width = 0.8)) +
  geom_linerange(aes(ymin = CI_lower, ymax = CI_upper), 
                 position = position_dodge(width = 0.8), na.rm = TRUE) +
  scale_x_continuous(breaks = seq.int(1, 10, 1)) +
  facet_wrap(~ estimate, scales = "free_y", ncol = 1) +
  labs(title = "Sensitivity analysis of GDM-birth weight mediation via isoform expression",
       subtitle = "Data from GUSTO (n=200)",
       x = "PCs (Principal Components)",
       y = "Estimate",
       color = "Annotation",
       fill = "Annotation") +
  theme_minimal() +
  theme(legend.position = "right") +
  scale_fill_manual(values = rev(c("seagreen", "#de77ae", "royalblue"))) + 
  scale_color_manual(values = rev(c("seagreen", "#de77ae", "royalblue")))
figs8a

ggsave("Figures/figs8a.pdf", figs8a, width = 7.5, height = 4.5, units = "in", dpi = 300)


################################################################################
# Figure S8D - PC mediation sensitivity analysis for GUSTO with GTEx tissues
################################################################################
PCsensitivity_d <- PCsensitivity %>%
  filter(Annotation %in% c("Placenta (lr-assembly)", "Adipose - Subcutaneous", "Cells - Cultured Fibroblasts"))

# Prepare data with estimate as category
PCsensitivity_d_long <- PCsensitivity_d %>%
  pivot_longer(cols = c(total_indirect, percent_variance),
               names_to = "estimate",
               values_to = "value") %>%
  mutate(estimate = factor(estimate, 
                           levels = c("total_indirect", "percent_variance"),
                           labels = c("Indirect Effect", "% Variance Explained")))

# Adjust CI columns for variance (set to NA since variance doesn't have CIs)
PCsensitivity_d_long <- PCsensitivity_d_long %>%
  mutate(CI_lower = ifelse(estimate == "Indirect Effect", CI_lower, NA_real_),
         CI_upper = ifelse(estimate == "Indirect Effect", CI_upper, NA_real_))

# Create data for conditional horizontal lines
hline_data <- data.frame(
  estimate = factor(c("Indirect Effect", "% Variance Explained"),
                    levels = c("Indirect Effect", "% Variance Explained")),
  yintercept = c(0, 75)
)

figs8d <- ggplot(PCsensitivity_d_long, aes(x = PCs, y = value, color = Annotation, fill = Annotation)) +
  geom_hline(data = hline_data, aes(yintercept = yintercept), linetype = "dotted") + 
  geom_point(size = 1.5, position = position_dodge(width = 0.8)) +
  geom_linerange(aes(ymin = CI_lower, ymax = CI_upper), 
                 position = position_dodge(width = 0.8), na.rm = TRUE) +
  scale_x_continuous(breaks = seq.int(1, 10, 1)) +
  facet_wrap(~ estimate, scales = "free_y", ncol = 1) +
  labs(title = "Sensitivity analysis of GDM-birth weight mediation via isoform expression",
       subtitle = "Data from GUSTO (n=200)",
       x = "PCs (Principal Components)",
       y = "Estimate",
       color = "Annotation",
       fill = "Annotation") +
  theme_minimal() +
  theme(legend.position = "right") +
  scale_fill_manual(values = rev(c("#ff68a1", "#e68613", "seagreen"))) + 
  scale_color_manual(values = rev(c("#ff68a1", "#e68613", "seagreen")))
figs8d

ggsave("Figures/figs8d.pdf", figs8d, width = 7.5, height = 4.5, units = "in", dpi = 300)




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
  
  results_indirect <- list()
  results_total <- list()
  
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
    
    # Extract total_indirect statistics
    ind_row <- param_est %>% filter(label == "total_indirect") %>%
      dplyr::select(est = est, ci.lower = ci.lower, ci.upper = ci.upper, pvalue = pvalue)
    
    # Extract total effect and proportion mediated statistics
    total_row <- param_est %>% filter(label == "total") %>%
      dplyr::select(est = est, ci.lower = ci.lower, ci.upper = ci.upper, pvalue = pvalue)
    
    prop_row <- param_est %>% filter(label == "prop_mediated") %>%
      dplyr::select(est = est, ci.lower = ci.lower, ci.upper = ci.upper, pvalue = pvalue)
    
    # Save indirect effect results
    results_indirect[[i]] <- data.frame(
      PCs = i,
      percent_variance = round(pca_cumvar[i] * 100, 2),
      total_indirect = ind_row$est,
      CI_lower = ind_row$ci.lower,
      CI_upper = ind_row$ci.upper,
      p = ind_row$pvalue
    )
    
    # Save total effect results
    results_total[[i]] <- data.frame(
      PCs = i,
      percent_variance = round(pca_cumvar[i] * 100, 2),
      total_effect = total_row$est,
      CI_lower = total_row$ci.lower,
      CI_upper = total_row$ci.upper,
      p = total_row$pvalue,
      prop_mediated = prop_row$est,
      prop_CI_lower = prop_row$ci.lower,
      prop_CI_upper = prop_row$ci.upper,
      prop_p = prop_row$pvalue
    )
  }
  
  final_results <- list(
    indirect = bind_rows(results_indirect),
    total = bind_rows(results_total)
  )
  
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
Gen3G.assembly.PCsensitivity.all <- PCmed.sensitivity.Gen3G(Gen3G.expr.assembly, cordat.Gen3G, 15)
Gen3G.assembly.PCsensitivity.totals <- Gen3G.assembly.PCsensitivity.all$total
Gen3G.assembly.PCsensitivity <- Gen3G.assembly.PCsensitivity.all$indirect

write.csv(Gen3G.assembly.PCsensitivity.totals,
          "Supplementary/PCsensitivity.totals.Gen3G.csv",row.names=F)
write.csv(Gen3G.assembly.PCsensitivity,
          "Supplementary/PCsensitivity.Gen3G.csv",row.names=F)

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
Gen3G.gencode.PCsensitivity.all <- PCmed.sensitivity.Gen3G(Gen3G.expr.gencode, cordat.Gen3G, 15)
Gen3G.gencode.PCsensitivity.totals <- Gen3G.gencode.PCsensitivity.all$total
Gen3G.gencode.PCsensitivity <- Gen3G.gencode.PCsensitivity.all$indirect

write.csv(Gen3G.gencode.PCsensitivity.totals,
          "Supplementary/PCsensitivity.totals.Gen3G.gencode.csv",row.names=F)
write.csv(Gen3G.gencode.PCsensitivity,
          "Supplementary/PCsensitivity.Gen3G.gencode.csv",row.names=F)

## Gen3G gencode_plus: PC-based mediation analysis
Gen3G.expr.gencode_plus <- data.frame(Gen3G.vst.gencode_plus)
Gen3G.expr.gencode_plus$transcript_id <- row.names(Gen3G.expr.gencode_plus)

# Filter to most significant transcripts for PC analysis
Gen3G.expr.gencode_plus <- Gen3G.expr.gencode_plus[Gen3G.expr.gencode_plus$transcript_id %in% 
                                                     results.gencode_plus.Gen3G.gdm[results.gencode_plus.Gen3G.gdm$pvalue < 0.00005, "transcript_id"], ]
Gen3G.expr.gencode_plus$transcript_id <- NULL

## Examine principal components 1-20
pca_result <- prcomp(t(as.matrix(Gen3G.expr.gencode_plus)))
variance_explained <- pca_result$sdev^2
percent_variance <- variance_explained / sum(variance_explained) * 100
percent_variance[1:20]
cumulative_variance <- cumsum(percent_variance)
cumulative_variance[1:20]

plot(percent_variance[1:20], type = "b", 
     xlab = "Principal Component", 
     ylab = "Percentage Variance Explained")

# Run PC sensitivity analysis for Gen3G gencode_plus
Gen3G.gencode_plus.PCsensitivity.all <- PCmed.sensitivity.Gen3G(Gen3G.expr.gencode_plus, cordat.Gen3G, 15)
Gen3G.gencode_plus.PCsensitivity.totals <- Gen3G.gencode_plus.PCsensitivity.all$total
Gen3G.gencode_plus.PCsensitivity <- Gen3G.gencode_plus.PCsensitivity.all$indirect

write.csv(Gen3G.gencode_plus.PCsensitivity.totals,
          "Supplementary/PCsensitivity.totals.Gen3G.gencode_plus.csv",row.names=F)
write.csv(Gen3G.gencode_plus.PCsensitivity,
          "Supplementary/PCsensitivity.Gen3G.gencode_plus.csv",row.names=F)

## Create sensitivity plot for Gen3G (Figure S8)
Gen3G.assembly.PCsensitivity$Annotation <- "lr-assembly"
Gen3G.gencode.PCsensitivity$Annotation <- "GENCODEv45"
Gen3G.gencode_plus.PCsensitivity$Annotation <- "GENCODE+"
PCsensitivity <- rbind(Gen3G.assembly.PCsensitivity,
                       Gen3G.gencode.PCsensitivity,
                       Gen3G.gencode_plus.PCsensitivity)
rownames(PCsensitivity) <- NULL

# Figure S8B - PC mediation sensitivity analysis for Gen3G
# Prepare data with estimate as category
PCsensitivity_long <- PCsensitivity %>%
  filter(PCs %in% 1:15) %>%
  pivot_longer(cols = c(total_indirect, percent_variance),
               names_to = "estimate",
               values_to = "value") %>%
  mutate(estimate = factor(estimate, 
                           levels = c("total_indirect", "percent_variance"),
                           labels = c("Indirect Effect", "% Variance Explained")))

# Adjust CI columns for variance (set to NA since variance doesn't have CIs)
PCsensitivity_long <- PCsensitivity_long %>%
  mutate(CI_lower = ifelse(estimate == "Indirect Effect", CI_lower, NA_real_),
         CI_upper = ifelse(estimate == "Indirect Effect", CI_upper, NA_real_))


################################################################################
# Figure S8B - PC mediation sensitivity analysis for Gen3G
################################################################################

# Create data for conditional horizontal lines
hline_data <- data.frame(
  estimate = factor(c("Indirect Effect", "% Variance Explained"),
                    levels = c("Indirect Effect", "% Variance Explained")),
  yintercept = c(0, 75)
)

figs8b <- ggplot(PCsensitivity_long, aes(x = PCs, y = value, color = Annotation, fill = Annotation)) +
  geom_hline(data = hline_data, aes(yintercept = yintercept), linetype = "dotted") + 
  geom_point(size = 1.5, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = CI_lower, ymax = CI_upper), 
                 position = position_dodge(width = 0.5), na.rm = TRUE) +
  scale_x_continuous(breaks = seq.int(1, 15, 1)) +
  facet_wrap(~ estimate, scales = "free_y", ncol = 1) +
  labs(title = "Sensitivity analysis of GDM-birth weight mediation via isoform expression",
       subtitle = "Data from Gen3G (n=152)",
       x = "PCs (Principal Components)",
       y = "Estimate",
       color = "Annotation",
       fill = "Annotation") +
  theme_minimal() +
  theme(legend.position = "right") +
  scale_fill_manual(values = c("royalblue", "#de77ae", "seagreen")) +
  scale_color_manual(values = c("royalblue", "#de77ae", "seagreen"))
figs8b

ggsave("Figures/figs8b.pdf", figs8b, width = 7.5, height = 4.5, units = "in", dpi = 300)



################################################################################
# Save results tables for running the remainder of analyses locally with internet access
################################################################################
save(
  list = c(
    "tx2g.gencode",
    "tx2g.assembly",
    "classif",
    "results.assembly.Gen3G.gdm",
    "results.gencode.Gen3G.gdm",
    "results.gencode_plus.Gen3G.gdm",
    "results.assembly.GUSTO.gdm",
    "results.gencode.GUSTO.gdm",
    "results.gencode_plus.GUSTO.gdm",
    "results.adipose.GUSTO.gdm",
    "results.fibroblasts.GUSTO.gdm",
    "results.assembly.Gen3G.gdm.gene",
    "results.gencode.Gen3G.gdm.gene",
    "results.gencode_plus.Gen3G.gdm.gene",
    "results.assembly.GUSTO.gdm.gene",
    "results.gencode.GUSTO.gdm.gene",
    "results.gencode_plus.GUSTO.gdm.gene",
    "results.adipose.GUSTO.gdm.gene",
    "results.fibroblasts.GUSTO.gdm.gene",
    "GUSTO.mediation.gdm",
    "Gen3G.mediation.gdm",
    "GUSTO.vst.assembly",
    "Gen3G.vst.assembly",
    "cordat.GUSTO",
    "cordat.Gen3G"
  ),
  file = "Supplementary/all_results_tables.RData"
)


################################################################################
# Fig3F - Direct & indirect effects across cohorts and annotations
################################################################################

load("Supplementary/all_results_tables.RData")

plot_mediation_effects <- function(pc_specs, cohort = c("GUSTO", "Gen3G")) {
  
  require(ggplot2)
  require(dplyr)
  
  # Validate cohort argument
  cohort <- match.arg(cohort, c("GUSTO", "Gen3G"), several.ok = TRUE)
  
  # pc_specs should be a named vector with elements for requested cohorts
  # e.g., c("GUSTO_GENCODEv45" = X, "GUSTO_GENCODE+" = X, "GUSTO_lr-assembly" = X)
  # or all 6 if both cohorts requested
  
  # Build expected order based on requested cohorts
  expected_order <- c()
  if ("GUSTO" %in% cohort) {
    expected_order <- c(expected_order, "GUSTO_GENCODEv45", "GUSTO_GENCODE+", "GUSTO_lr-assembly")
  }
  if ("Gen3G" %in% cohort) {
    expected_order <- c(expected_order, "Gen3G_GENCODEv45", "Gen3G_GENCODE+", "Gen3G_lr-assembly")
  }
  
  if (!all(names(pc_specs) == expected_order)) {
    stop("pc_specs must have names in this exact order: ", paste(expected_order, collapse = ", "))
  }
  
  # Read in requested cohort files
  datasets <- list()
  
  if ("GUSTO" %in% cohort) {
    GUSTO.assembly.indirect <- read.csv("Supplementary/PCsensitivity.GUSTO.csv")
    GUSTO.assembly.totals <- read.csv("Supplementary/PCsensitivity.totals.GUSTO.csv")
    GUSTO.gencode.indirect <- read.csv("Supplementary/PCsensitivity.GUSTO.gencode.csv")
    GUSTO.gencode.totals <- read.csv("Supplementary/PCsensitivity.totals.GUSTO.gencode.csv")
    GUSTO.gencode_plus.indirect <- read.csv("Supplementary/PCsensitivity.GUSTO.gencode_plus.csv")
    GUSTO.gencode_plus.totals <- read.csv("Supplementary/PCsensitivity.totals.GUSTO.gencode_plus.csv")
    
    GUSTO.gencode.indirect$Cohort <- "GUSTO"
    GUSTO.gencode.indirect$Annotation <- "GENCODEv45"
    GUSTO.gencode_plus.indirect$Cohort <- "GUSTO"
    GUSTO.gencode_plus.indirect$Annotation <- "GENCODE+"
    GUSTO.assembly.indirect$Cohort <- "GUSTO"
    GUSTO.assembly.indirect$Annotation <- "lr-assembly"
    
    GUSTO.gencode.totals$Cohort <- "GUSTO"
    GUSTO.gencode.totals$Annotation <- "GENCODEv45"
    GUSTO.gencode_plus.totals$Cohort <- "GUSTO"
    GUSTO.gencode_plus.totals$Annotation <- "GENCODE+"
    GUSTO.assembly.totals$Cohort <- "GUSTO"
    GUSTO.assembly.totals$Annotation <- "lr-assembly"
    
    GUSTO.gencode.indirect <- GUSTO.gencode.indirect %>% filter(PCs == pc_specs["GUSTO_GENCODEv45"])
    GUSTO.gencode_plus.indirect <- GUSTO.gencode_plus.indirect %>% filter(PCs == pc_specs["GUSTO_GENCODE+"])
    GUSTO.assembly.indirect <- GUSTO.assembly.indirect %>% filter(PCs == pc_specs["GUSTO_lr-assembly"])
    GUSTO.gencode.totals <- GUSTO.gencode.totals %>% filter(PCs == pc_specs["GUSTO_GENCODEv45"])
    GUSTO.gencode_plus.totals <- GUSTO.gencode_plus.totals %>% filter(PCs == pc_specs["GUSTO_GENCODE+"])
    GUSTO.assembly.totals <- GUSTO.assembly.totals %>% filter(PCs == pc_specs["GUSTO_lr-assembly"])
    
    datasets$GUSTO <- list(
      indirect = bind_rows(GUSTO.gencode.indirect, GUSTO.gencode_plus.indirect, GUSTO.assembly.indirect),
      totals = bind_rows(GUSTO.gencode.totals, GUSTO.gencode_plus.totals, GUSTO.assembly.totals)
    )
  }
  
  if ("Gen3G" %in% cohort) {
    Gen3G.assembly.indirect <- read.csv("Supplementary/PCsensitivity.Gen3G.csv")
    Gen3G.assembly.totals <- read.csv("Supplementary/PCsensitivity.totals.Gen3G.csv")
    Gen3G.gencode.indirect <- read.csv("Supplementary/PCsensitivity.Gen3G.gencode.csv")
    Gen3G.gencode.totals <- read.csv("Supplementary/PCsensitivity.totals.Gen3G.gencode.csv")
    Gen3G.gencode_plus.indirect <- read.csv("Supplementary/PCsensitivity.Gen3G.gencode_plus.csv")
    Gen3G.gencode_plus.totals <- read.csv("Supplementary/PCsensitivity.totals.Gen3G.gencode_plus.csv")
    
    Gen3G.gencode.indirect$Cohort <- "Gen3G"
    Gen3G.gencode.indirect$Annotation <- "GENCODEv45"
    Gen3G.gencode_plus.indirect$Cohort <- "Gen3G"
    Gen3G.gencode_plus.indirect$Annotation <- "GENCODE+"
    Gen3G.assembly.indirect$Cohort <- "Gen3G"
    Gen3G.assembly.indirect$Annotation <- "lr-assembly"
    
    Gen3G.gencode.totals$Cohort <- "Gen3G"
    Gen3G.gencode.totals$Annotation <- "GENCODEv45"
    Gen3G.gencode_plus.totals$Cohort <- "Gen3G"
    Gen3G.gencode_plus.totals$Annotation <- "GENCODE+"
    Gen3G.assembly.totals$Cohort <- "Gen3G"
    Gen3G.assembly.totals$Annotation <- "lr-assembly"
    
    Gen3G.gencode.indirect <- Gen3G.gencode.indirect %>% filter(PCs == pc_specs["Gen3G_GENCODEv45"])
    Gen3G.gencode_plus.indirect <- Gen3G.gencode_plus.indirect %>% filter(PCs == pc_specs["Gen3G_GENCODE+"])
    Gen3G.assembly.indirect <- Gen3G.assembly.indirect %>% filter(PCs == pc_specs["Gen3G_lr-assembly"])
    Gen3G.gencode.totals <- Gen3G.gencode.totals %>% filter(PCs == pc_specs["Gen3G_GENCODEv45"])
    Gen3G.gencode_plus.totals <- Gen3G.gencode_plus.totals %>% filter(PCs == pc_specs["Gen3G_GENCODE+"])
    Gen3G.assembly.totals <- Gen3G.assembly.totals %>% filter(PCs == pc_specs["Gen3G_lr-assembly"])
    
    datasets$Gen3G <- list(
      indirect = bind_rows(Gen3G.gencode.indirect, Gen3G.gencode_plus.indirect, Gen3G.assembly.indirect),
      totals = bind_rows(Gen3G.gencode.totals, Gen3G.gencode_plus.totals, Gen3G.assembly.totals)
    )
  }
  
  # Combine all indirect and total effects
  all_indirect <- bind_rows(lapply(datasets, function(x) x$indirect))
  all_totals <- bind_rows(lapply(datasets, function(x) x$totals))
  
  # Cap proportion mediated CIs at -1 and 1
  all_totals <- all_totals %>%
    dplyr::mutate(
      prop_CI_lower = pmax(prop_CI_lower, -1),
      prop_CI_upper = pmin(prop_CI_upper, 1)
    )
  
  # Prepare data for plotting
  # Indirect effects
  indirect_plot <- all_indirect %>%
    dplyr::select(Cohort, Annotation, PCs, estimate = total_indirect, CI_lower, CI_upper, percent_variance) %>%
    dplyr::mutate(Effect = "Indirect")
  
  # Total effects - keep only one per cohort (they're the same across annotations)
  total_plot <- all_totals %>%
    dplyr::group_by(Cohort) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(Cohort, Annotation, PCs, estimate = total_effect, CI_lower, CI_upper, percent_variance) %>%
    dplyr::mutate(Effect = "Total",
                  Annotation = "Total")  # Special annotation label
  
  # Proportion mediated
  prop_plot <- all_totals %>%
    dplyr::select(Cohort, Annotation, PCs, estimate = prop_mediated, 
                  CI_lower = prop_CI_lower, CI_upper = prop_CI_upper, percent_variance) %>%
    dplyr::mutate(Effect = "Proportion Mediated")
  
  # Combine for plotting
  plot_data <- bind_rows(total_plot, indirect_plot, prop_plot)
  
  # Set factor levels for ordering
  plot_data$Annotation <- factor(plot_data$Annotation, 
                                 levels = c("GENCODEv45", "GENCODE+", "lr-assembly", "Total"))
  plot_data$Cohort <- factor(plot_data$Cohort, levels = c("GUSTO", "Gen3G"))
  plot_data$Effect <- factor(plot_data$Effect, levels = rev(c("Total", "Indirect", "Proportion Mediated")))
  
  # Create color scale - black for Total, colors for others
  color_values <- c("GENCODEv45" = "royalblue", 
                    "GENCODE+" = "#de77ae", 
                    "lr-assembly" = "seagreen",
                    "Total" = "black")
  
  plot_data[plot_data$Effect=="Proportion Mediated"&plot_data$Cohort=="GUSTO","CI_lower"] <- plot_data[plot_data$Effect=="Proportion Mediated"&plot_data$Cohort=="GUSTO","CI_lower"]+.01
  
  if(length(cohort)>1){cohort.lab="all cohorts"}else{
    cohort.lab <- cohort
  }
  if(cohort.lab=="GUSTO"){cohort.lab <- "GUSTO (n=200)"}
  if(cohort.lab=="Gen3G"){cohort.lab <- "Gen3G (n=152)"}
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = estimate, y = Effect, color = Annotation, group = Annotation)) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_point(size = 2, position = position_dodge(width = 0.6)) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), 
                   position = position_dodge(width = 0.6), height = 0.2, linewidth = 0.5) +
    scale_color_manual(values = color_values,
                       breaks = c("GENCODEv45", "GENCODE+", "lr-assembly")) +  # Exclude "Total" from legend
    labs(title = "GDM-birth weight mediation via placental isoform expression PCs",
         subtitle = bquote("Selected PCs explain" >= "75% of variance; data from " * .(paste(cohort.lab, collapse = " and "))),
         x = "Effect Estimate",
         y = "",
         color = "Annotation") +
    theme_minimal() +
    theme(legend.position = "bottom",
          panel.grid.major.y = element_blank(),
          strip.background = element_rect(fill = "white", color = "black"),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(hjust = 0),
          plot.subtitle = element_text(hjust = 0),
          plot.title.position = "plot",
          plot.caption.position = "plot")
  
  # Add faceting only if multiple cohorts
  if (length(cohort) > 1) {
    p <- p + facet_wrap(~Cohort, nrow = 2)
  }
  
  return(p)
}

# pc_specs <- c("GUSTO_GENCODEv45" = 8, "GUSTO_GENCODE+" = 6, "GUSTO_lr-assembly" = 5,
#               "Gen3G_GENCODEv45" = 14, "Gen3G_GENCODE+" = 7, "Gen3G_lr-assembly" = 6)

pc_specs <- c("GUSTO_GENCODEv45" = 8, "GUSTO_GENCODE+" = 6, "GUSTO_lr-assembly" = 5)

fig3f <- plot_mediation_effects(pc_specs,"GUSTO")
fig3f

ggsave("Figures/fig3f.pdf", fig3f,  width = 5.76, height = 3.68, units = "in", dpi = 300)




################################################################################
# Gene Ontology Enrichment Analysis Setup
################################################################################

# Create directories for enrichment analysis results
dir.create("EnrichR_Results/EnrichR_Results_DTE", showWarnings = FALSE)
dir.create("EnrichR_Results/EnrichR_Results_DTE/Plots", showWarnings = FALSE)
dir.create("EnrichR_Results/EnrichR_Results_DGE", showWarnings = FALSE)
dir.create("EnrichR_Results/EnrichR_Results_DGE/Plots", showWarnings = FALSE)

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
  
  # Create the lollipop plot — no axis break
  p <- ggplot(df, aes(x = neglog10p, y = Ontology)) +
    geom_segment(aes(x = 0, xend = neglog10p,
                     y = Ontology, yend = Ontology),
                 color = "gray50", linewidth = 0.7) +
    geom_point(aes(size = Odds.Ratio, color = Database),
               alpha = 0.8) +
    scale_size_continuous(range = c(3, 8),
                          name = "Odds Ratio") +
    scale_color_manual(values = db_palette) +
    labs(
      title = paste("Top 15 Enriched Terms -", dataset_name),
      subtitle = "Point size = Odds Ratio, X-axis = -log10(FDR)",
      x = "-log10(Adjusted P-value)",
      y = ""
    ) +
    theme_bw() +
    theme(
      plot.title    = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10),
      axis.text.y   = element_text(size = 8),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.position = "top",
      legend.title    = element_blank()
    )
  
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
  
  # Use continuous x scale as usual
  p <- p + scale_x_continuous(
    expand = expansion(mult = 0),
    limits = c(0, max(df$neglog10p + 1))
  )
  
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
                                      "EnrichR_Results/EnrichR_Results_top50isorichgenes")
top50_summary <- create_summary_table(top50_enrichr,
                                      NULL,
                                      "EnrichR_Results/EnrichR_Results_top50isorichgenes/results.csv")

# Extract molecular function results for top 50 isoform-rich genes
top50MF <- top50_summary[top50_summary$Database == "GO_Molecular_Function_2023", ]
write.csv(top50MF, "top50isorichgenes_GOMF.csv", row.names = FALSE)

################################################################################
# Differential Transcript Expression (DTE) Enrichment Analysis
# Includes Figures S7 C-D
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
                                                       "EnrichR_Results/EnrichR_Results_DTE")

enrichr_gencode_results_Gen3G <- run_enrichr_analysis(symbols_sig_gencode_Gen3G, 
                                                      symbols_universe_gencode, 
                                                      "Gen3G_gencode_sig", dbs,
                                                      "EnrichR_Results/EnrichR_Results_DTE")

enrichr_assembly_results_GUSTO <- run_enrichr_analysis(symbols_sig_assembly_GUSTO, 
                                                       symbols_universe_assembly, 
                                                       "GUSTO_assembly_sig", dbs,
                                                       "EnrichR_Results/EnrichR_Results_DTE")

enrichr_gencode_results_GUSTO <- run_enrichr_analysis(symbols_sig_gencode_GUSTO, 
                                                      symbols_universe_gencode, 
                                                      "GUSTO_gencode_sig", dbs,
                                                      "EnrichR_Results/EnrichR_Results_DTE")

## Create combined summary tables for DTE
all_results_Gen3G <- create_summary_table(enrichr_assembly_results_Gen3G,
                                          enrichr_gencode_results_Gen3G,
                                          "EnrichR_Results/EnrichR_Results_DTE/results_Gen3G_DTE.csv")
all_results_Gen3G.DTE <- all_results_Gen3G

all_results_GUSTO <- create_summary_table(enrichr_assembly_results_GUSTO,
                                          enrichr_gencode_results_GUSTO,
                                          "EnrichR_Results/EnrichR_Results_DTE/results_GUSTO_DTE.csv")
all_results_GUSTO.DTE <- all_results_GUSTO

## Create lollipop plots for DTE visualization
# Get list of datasets
datasets_Gen3G <- unique(all_results_Gen3G$Dataset)
datasets_GUSTO <- unique(all_results_GUSTO$Dataset)

# Create lollipop plots for each dataset
make_lolipops_perdataset(datasets_Gen3G, all_results_Gen3G,
                         "EnrichR_Results/EnrichR_Results_DTE/Plots/",
                         "_lollipop_DTE.pdf")

make_lolipops_perdataset(datasets_GUSTO, all_results_GUSTO,
                         "EnrichR_Results/EnrichR_Results_DTE/Plots/",
                         "_lollipop_DTE.pdf")

################################################################################
# Differential Gene Expression (DGE) Enrichment Analysis
################################################################################

## Prepare gene-level results
# Subset columns and set Confirmation.P as padj
results.assembly.Gen3G.gdm.gene$gene_id <- sapply(strsplit(results.assembly.Gen3G.gdm.gene$gene_id, '[.]'), function(x) x[1])

results.gencode.Gen3G.gdm.gene$gene_id <- sapply(strsplit(results.gencode.Gen3G.gdm.gene$gene_id, '[.]'), function(x) x[1])

results.assembly.GUSTO.gdm.gene$gene_id <- sapply(strsplit(results.assembly.GUSTO.gdm.gene$gene_id, '[.]'), function(x) x[1])

results.gencode.GUSTO.gdm.gene$gene_id <- sapply(strsplit(results.gencode.GUSTO.gdm.gene$gene_id, '[.]'), function(x) x[1])

## Identify genes with padj < 0.1
sig_transcripts_assembly_Gen3G <- results.assembly.Gen3G.gdm.gene[results.assembly.Gen3G.gdm.gene$padj < 0.1, ]
sig_transcripts_gencode_Gen3G <- results.gencode.Gen3G.gdm.gene[results.gencode.Gen3G.gdm.gene$padj < 0.1, ]
sig_transcripts_assembly_GUSTO <- results.assembly.GUSTO.gdm.gene[results.assembly.GUSTO.gdm.gene$padj < 0.1, ]
sig_transcripts_gencode_GUSTO <- results.gencode.GUSTO.gdm.gene[results.gencode.GUSTO.gdm.gene$padj < 0.1, ]

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
                                                       "EnrichR_Results/EnrichR_Results_DGE")

enrichr_gencode_results_Gen3G <- run_enrichr_analysis(symbols_sig_gencode_Gen3G, 
                                                      symbols_universe_gencode, 
                                                      "Gen3G_gencode_sig", dbs,
                                                      "EnrichR_Results/EnrichR_Results_DGE")

enrichr_assembly_results_GUSTO <- run_enrichr_analysis(symbols_sig_assembly_GUSTO, 
                                                       symbols_universe_assembly, 
                                                       "GUSTO_assembly_sig", dbs,
                                                       "EnrichR_Results/EnrichR_Results_DGE")

enrichr_gencode_results_GUSTO <- run_enrichr_analysis(symbols_sig_gencode_GUSTO, 
                                                      symbols_universe_gencode, 
                                                      "GUSTO_gencode_sig", dbs,
                                                      "EnrichR_Results/EnrichR_Results_DGE")

## Create combined summary tables for DGE
all_results_Gen3G <- create_summary_table(enrichr_assembly_results_Gen3G,
                                          enrichr_gencode_results_Gen3G,
                                          "EnrichR_Results/EnrichR_Results_DGE/results_Gen3G_DGE.csv")
all_results_Gen3G.DGE <- all_results_Gen3G

all_results_GUSTO <- create_summary_table(enrichr_assembly_results_GUSTO,
                                          enrichr_gencode_results_GUSTO,
                                          "EnrichR_Results/EnrichR_Results_DGE/results_GUSTO_DGE.csv")
all_results_GUSTO.DGE <- all_results_GUSTO

## Create lollipop plots for DGE visualization
# Get list of datasets
datasets_Gen3G <- unique(all_results_Gen3G$Dataset)
datasets_GUSTO <- unique(all_results_GUSTO$Dataset)

# Create lollipop plots for each dataset
make_lolipops_perdataset(datasets_Gen3G, all_results_Gen3G,
                         "EnrichR_Results/EnrichR_Results_DGE/Plots/",
                         "_lollipop_DGE.pdf")

make_lolipops_perdataset(datasets_GUSTO, all_results_GUSTO,
                         "EnrichR_Results/EnrichR_Results_DGE/Plots/",
                         "_lollipop_DGE.pdf")




################################################################################
# Mediation Results Visualization and Functional Annotation
################################################################################

##### ----- Run locally ----- #####

# Set up biomart for gene annotation
## mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

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
# GUSTO Mediation Results Annotation
################################################################################

## Process GUSTO significant mediation results
GUSTO.mediation.gdm.sig <- GUSTO.mediation.gdm[GUSTO.mediation.gdm$pvalue < 0.05, ]
GUSTO.mediation.gdm.sig <- GUSTO.mediation.gdm.sig[GUSTO.mediation.gdm.sig$annotation == "lr-assembly", ]
GUSTO.mediation.gdm.sig <- GUSTO.mediation.gdm.sig[order(GUSTO.mediation.gdm.sig$estimate_abs, decreasing = TRUE), ]
GUSTO.mediation.gdm.sig$geneID <- sapply(strsplit(GUSTO.mediation.gdm.sig$gene_id, '[.]'), function(x) x[1])

# Map to Entrez IDs for GO analysis
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = GUSTO.mediation.gdm.sig$geneID,
                     column = "ENTREZID",
                     keytype = "ENSEMBL",
                     multiVals = "first")
names(entrez_ids) <- GUSTO.mediation.gdm.sig$geneID

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
names(GUSTO.mediation.gdm.sig)[19] <- "Function"
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

GUSTO.mediation.gdm.sig <- left_join(GUSTO.mediation.gdm.sig, gene_symbols)
GUSTO.mediation.gdm.sig$label <- paste0(GUSTO.mediation.gdm.sig$symbol, "\n(",
                                        GUSTO.mediation.gdm.sig$transcript_id, ")")

# Prepare data for visualization
GUSTO.mediation.gdm.sig <- GUSTO.mediation.gdm.sig[order(GUSTO.mediation.gdm.sig$estimate_abs), ]

GUSTO.mediation.gdm.sig$Function <- factor(GUSTO.mediation.gdm.sig$Function,
                                           levels = c("Cell Organization", "Metabolism",
                                                      "Development", "Immune", "Regulation",
                                                      "Unknown"))

GUSTO.mediation.gdm.sig <- GUSTO.mediation.gdm.sig[order(GUSTO.mediation.gdm.sig$estimate_abs, decreasing = TRUE), ]
GUSTO.mediation.gdm.sig.a <- GUSTO.mediation.gdm.sig


################################################################################
# Gen3G Mediation Results Annotation
################################################################################

## Process Gen3G significant mediation results
Gen3G.mediation.gdm.sig <- Gen3G.mediation.gdm[Gen3G.mediation.gdm$pvalue < 0.05, ]
Gen3G.mediation.gdm.sig <- Gen3G.mediation.gdm.sig[Gen3G.mediation.gdm.sig$annotation == "lr-assembly", ]
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
names(Gen3G.mediation.gdm.sig)[19] <- "Function"
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
# gene_symbols[5,2] <- "Lnc-ISOC1-1"

Gen3G.mediation.gdm.sig <- left_join(Gen3G.mediation.gdm.sig, gene_symbols)
Gen3G.mediation.gdm.sig$label <- paste0(Gen3G.mediation.gdm.sig$symbol, "\n(",
                                        Gen3G.mediation.gdm.sig$transcript_id, ")")

# Prepare data for visualization
Gen3G.mediation.gdm.sig <- Gen3G.mediation.gdm.sig[order(Gen3G.mediation.gdm.sig$estimate_abs), ]

Gen3G.mediation.gdm.sig$Function <- factor(Gen3G.mediation.gdm.sig$Function,
                                           levels = c("Cell Organization", "Metabolism",
                                                      "Development", "Immune", "Regulation",
                                                      "Unknown"))

Gen3G.mediation.gdm.sig <- Gen3G.mediation.gdm.sig[!duplicated(Gen3G.mediation.gdm.sig),]

# Get top 10 lr-assembly transcripts from Gen3G
Gen3G.mediation.gdm.sig <- Gen3G.mediation.gdm.sig[order(Gen3G.mediation.gdm.sig$estimate_abs, decreasing = TRUE), ]
Gen3G.mediation.gdm.sig <- Gen3G.mediation.gdm.sig[1:10, ]
Gen3G.mediation.gdm.sig <- Gen3G.mediation.gdm.sig[!is.na(Gen3G.mediation.gdm.sig$transcript_id), ]
Gen3G.mediation.gdm.sig.a <- Gen3G.mediation.gdm.sig


################################################################################
# Create union of transcripts and filter both datasets
################################################################################

# Get union of transcript IDs
union_transcripts <- union(GUSTO.mediation.gdm.sig.a$transcript_id, 
                           Gen3G.mediation.gdm.sig.a$transcript_id)

##### ------------------ Run on HPC --------------------------------- ####

union_transcripts <- c(
  "ENST00000368130.9", "ESPRESSO:chr12:7288:0", "ESPRESSO:chr19:15403:3",
  "ENST00000351205.8", "ENST00000463456.5", "ENST00000455635.1",
  "ENST00000502610.1", "ENST00000595930.1", "ESPRESSO:chr17:12985:70",
  "ESPRESSO:chr11:6026:0", "ENST00000519548.5", "ENST00000460021.1",
  "ENST00000650735.1", "ESPRESSO:chr17:12985:62", "ESPRESSO:chr20:18555:1",
  "ENST00000253122.10", "ENST00000239223.4", "ENST00000369583.4",
  "ENST00000368285.8", "ENST00000466677.1", "ENST00000664340.1",
  "ENST00000544052.6", "ENST00000258301.6"
)

## GUSTO lr-assembly: Mediation analysis for union transcripts
GUSTO.expr.assembly.union <- data.frame(GUSTO.vst.assembly)
GUSTO.expr.assembly.union$transcript_id <- row.names(GUSTO.expr.assembly.union)

# Filter to union transcripts only
GUSTO.expr.assembly.union <- GUSTO.expr.assembly.union[GUSTO.expr.assembly.union$transcript_id %in% union_transcripts, ]

# Reshape data for mediation analysis
GUSTO.expr.assembly.union <- GUSTO.expr.assembly.union %>%
  pivot_longer(
    cols = -transcript_id,
    names_to = "SampleID",
    values_to = "TXE")

# Add covariate data
GUSTO.expr.assembly.union <- left_join(GUSTO.expr.assembly.union, 
                                       cordat.GUSTO[, c("SampleID", "gdm", "birth_weight", "sex", "GA")])
txids <- unique(GUSTO.expr.assembly.union$transcript_id)

# Parallel mediation analysis for each transcript
GUSTO.mediation.gdm.sig.a.union <- foreach(tx = 1:length(txids), 
                                           .packages = c("mediation"),
                                           .combine = rbind) %dopar% {
                                             
                                             data <- GUSTO.expr.assembly.union[GUSTO.expr.assembly.union$transcript_id == txids[tx], ]
                                             mediator_model <- lm(TXE ~ gdm + sex + GA, data = data)
                                             outcome_model <- lm(birth_weight ~ gdm + TXE + sex + GA, data = data)
                                             
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

row.names(GUSTO.mediation.gdm.sig.a.union) <- NULL
GUSTO.mediation.gdm.sig.a.union$estimate_abs <- abs(GUSTO.mediation.gdm.sig.a.union$estimate)
GUSTO.mediation.gdm.sig.a.union <- left_join(GUSTO.mediation.gdm.sig.a.union, tx2g.assembly)
GUSTO.mediation.gdm.sig.a.union$annotation <- "lr-assembly"

## Gen3G lr-assembly: Mediation analysis for union transcripts
Gen3G.expr.assembly.union <- data.frame(Gen3G.vst.assembly)
Gen3G.expr.assembly.union$transcript_id <- row.names(Gen3G.expr.assembly.union)

# Filter to union transcripts only
Gen3G.expr.assembly.union <- Gen3G.expr.assembly.union[Gen3G.expr.assembly.union$transcript_id %in% union_transcripts, ]

# Reshape data for mediation analysis
Gen3G.expr.assembly.union <- Gen3G.expr.assembly.union %>%
  pivot_longer(
    cols = -transcript_id,
    names_to = "Run",
    values_to = "TXE")

# Add covariate data
Gen3G.expr.assembly.union <- left_join(Gen3G.expr.assembly.union, 
                                       cordat.Gen3G[, c("Run", "gdm", "birth_weight", "GA")])
txids <- unique(Gen3G.expr.assembly.union$transcript_id)

# Parallel mediation analysis for each transcript
Gen3G.mediation.gdm.sig.a.union <- foreach(tx = 1:length(txids), 
                                           .packages = c("mediation"),
                                           .combine = rbind) %dopar% {
                                             
                                             data <- Gen3G.expr.assembly.union[Gen3G.expr.assembly.union$transcript_id == txids[tx], ]
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

row.names(Gen3G.mediation.gdm.sig.a.union) <- NULL
Gen3G.mediation.gdm.sig.a.union$estimate_abs <- abs(Gen3G.mediation.gdm.sig.a.union$estimate)
Gen3G.mediation.gdm.sig.a.union <- left_join(Gen3G.mediation.gdm.sig.a.union, tx2g.assembly, by = "transcript_id")
Gen3G.mediation.gdm.sig.a.union$annotation <- "lr-assembly"

save(list=c("Gen3G.mediation.gdm.sig.a.union","GUSTO.mediation.gdm.sig.a.union"),
     file="Supplementary/mediation_union_transcripts.RData")


##### ------------------ Return to local machine --------------------------------- ####

load("Supplementary/mediation_union_transcripts.RData")

union_mapping <- rbind(GUSTO.mediation.gdm.sig.a,Gen3G.mediation.gdm.sig.a)
union_mapping <- union_mapping[,c("transcript_id","symbol","label","Function")]
union_mapping <- union_mapping[!duplicated(union_mapping),]

# Get info for missing transcripts
missing <- data.frame(transcript_id=setdiff(GUSTO.mediation.gdm.sig.a.union$transcript_id,union_mapping$transcript_id))
missing <- left_join(missing,tx2g.assembly)

missing$geneID <- sapply(strsplit(missing$gene_id, '[.]'), function(x) x[1])

entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = missing$geneID,
                     column = "ENTREZID",
                     keytype = "ENSEMBL",
                     multiVals = "first")

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

go_bp <- go_bp %>%
  mutate(priority_rank = match(top, priority)) %>%
  group_by(geneID) %>%
  slice_min(order_by = priority_rank, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(geneID, top)

missing <- left_join(missing, go_bp)
names(missing)[ncol(missing)] <- "Function"
missing[is.na(missing$Function), "Function"] <- "Unknown"

gene_symbols <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = missing$geneID,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

gene_symbols <- data.frame(geneID = missing$geneID, symbol = gene_symbols)
missing <- left_join(missing, gene_symbols)
missing$label <- paste0(missing$symbol, "\n(", missing$transcript_id, ")")

union_mapping <- rbind(union_mapping,missing[,c("transcript_id","symbol","label","Function")])

# Top 10 per cohort
keep_ids <- c(
  "ESPRESSO:chr19:15403:3",
  "ENST00000466677.1",
  "ENST00000258301.6",
  "ENST00000351205.8",
  "ENST00000253122.10",
  "ENST00000368285.8",
  "ENST00000544052.6",
  "ENST00000595930.1",
  "ENST00000502610.1",
  "ENST00000463456.5",
  "ENST00000664340.1",
  "ENST00000455635.1",
  "ESPRESSO:chr20:18555:1",
  "ENST00000369583.4",
  "ENST00000239223.4",
  "ESPRESSO:chr17:12985:70",
  "ESPRESSO:chr17:12985:62",
  "ESPRESSO:chr11:6026:0",
  "ESPRESSO:chr12:7288:0",
  "ENST00000368130.9"
)

union_mapping <- union_mapping[union_mapping$transcript_id%in%keep_ids,]



# Check for missing data and update

keep_cols <- names(GUSTO.mediation.gdm.sig.a.union)

### GUSTO
GUSTO.mediation.gdm.sig.a.union <- GUSTO.mediation.gdm.sig.a.union[!GUSTO.mediation.gdm.sig.a.union$transcript_id%in%GUSTO.mediation.gdm.sig.a$transcript_id,]
GUSTO.union <- union(GUSTO.mediation.gdm.sig.a$transcript_id,GUSTO.mediation.gdm.sig.a.union$transcript_id)
GUSTO.missing <- setdiff(keep_ids,GUSTO.union)

GUSTO.expr.assembly.missing <- data.frame(GUSTO.vst.assembly)
GUSTO.expr.assembly.missing$transcript_id <- row.names(GUSTO.expr.assembly.missing)

# Filter to GUSTO.missing transcripts only
GUSTO.expr.assembly.missing <- GUSTO.expr.assembly.missing[
  GUSTO.expr.assembly.missing$transcript_id %in% GUSTO.missing, ]

# Reshape data for mediation analysis
GUSTO.expr.assembly.missing <- GUSTO.expr.assembly.missing %>%
  pivot_longer(
    cols = -transcript_id,
    names_to = "SampleID",
    values_to = "TXE")

# Add covariate data
GUSTO.expr.assembly.missing <- left_join(
  GUSTO.expr.assembly.missing,
  cordat.GUSTO[, c("SampleID", "gdm", "birth_weight", "sex", "GA")]
)

txids <- unique(GUSTO.expr.assembly.missing$transcript_id)

# Parallel mediation analysis for each transcript
GUSTO.mediation.gdm.sig.a.missing <- foreach(
  tx = 1:length(txids),
  .packages = c("mediation"),
  .combine = rbind
) %dopar% {
  
  data <- GUSTO.expr.assembly.missing[
    GUSTO.expr.assembly.missing$transcript_id == txids[tx], ]
  
  mediator_model <- lm(TXE ~ gdm + sex + GA, data = data)
  outcome_model  <- lm(birth_weight ~ gdm + TXE + sex + GA, data = data)
  
  tryCatch({
    mediate_result <- mediate(
      mediator_model,
      outcome_model,
      treat = "gdm",
      mediator = "TXE",
      boot = TRUE,
      sims = 100
    )
    
    data.frame(
      transcript_id = txids[tx],
      estimate      = mediate_result$d0,
      lower_95CI    = mediate_result$d0.ci[1],
      upper_95CI    = mediate_result$d0.ci[2],
      pvalue        = mediate_result$d0.p
    )
    
  }, error = function(e) {
    
    warning(paste("Error in transcript", txids[tx], ":", e$message))
    
    data.frame(
      transcript_id = txids[tx],
      estimate      = NA,
      lower_95CI    = NA,
      upper_95CI    = NA,
      pvalue        = NA
    )
  })
}

row.names(GUSTO.mediation.gdm.sig.a.missing) <- NULL
GUSTO.mediation.gdm.sig.a.missing$estimate_abs <- 
  abs(GUSTO.mediation.gdm.sig.a.missing$estimate)

GUSTO.mediation.gdm.sig.a.missing <- 
  left_join(GUSTO.mediation.gdm.sig.a.missing, tx2g.assembly)

GUSTO.mediation.gdm.sig.a.missing$annotation <- "lr-assembly"

GUSTO.mediation.gdm.sig.a.union <- rbind(GUSTO.mediation.gdm.sig.a[,keep_cols],
                                         GUSTO.mediation.gdm.sig.a.union[,keep_cols],
                                         GUSTO.mediation.gdm.sig.a.missing[,keep_cols])
GUSTO.mediation.gdm.sig.a.union <- GUSTO.mediation.gdm.sig.a.union[GUSTO.mediation.gdm.sig.a.union$transcript_id%in%keep_ids,]



### Gen3G
Gen3G.mediation.gdm.sig.a.union <- Gen3G.mediation.gdm.sig.a.union[!Gen3G.mediation.gdm.sig.a.union$transcript_id%in%Gen3G.mediation.gdm.sig.a$transcript_id,]
Gen3G.union <- union(Gen3G.mediation.gdm.sig.a$transcript_id,Gen3G.mediation.gdm.sig.a.union$transcript_id)
Gen3G.missing <- setdiff(keep_ids,Gen3G.union)

Gen3G.expr.assembly.missing <- data.frame(Gen3G.vst.assembly)
Gen3G.expr.assembly.missing$transcript_id <- row.names(Gen3G.expr.assembly.missing)

# Filter to Gen3G.missing transcripts only
Gen3G.expr.assembly.missing <- Gen3G.expr.assembly.missing[
  Gen3G.expr.assembly.missing$transcript_id %in% Gen3G.missing, ]

# Reshape data for mediation analysis
Gen3G.expr.assembly.missing <- Gen3G.expr.assembly.missing %>%
  pivot_longer(
    cols = -transcript_id,
    names_to = "SampleID",
    values_to = "TXE")

# Add covariate data
names(Gen3G.expr.assembly.missing)[2] <- "Run"
Gen3G.expr.assembly.missing <- left_join(
  Gen3G.expr.assembly.missing,
  cordat.Gen3G[, c("Run", "gdm", "birth_weight", "sex", "GA")]
)

txids <- unique(Gen3G.expr.assembly.missing$transcript_id)

# Parallel mediation analysis for each transcript
Gen3G.mediation.gdm.sig.a.missing <- foreach(
  tx = 1:length(txids),
  .packages = c("mediation"),
  .combine = rbind
) %dopar% {
  
  data <- Gen3G.expr.assembly.missing[
    Gen3G.expr.assembly.missing$transcript_id == txids[tx], ]
  
  mediator_model <- lm(TXE ~ gdm + sex + GA, data = data)
  outcome_model  <- lm(birth_weight ~ gdm + TXE + sex + GA, data = data)
  
  tryCatch({
    mediate_result <- mediate(
      mediator_model,
      outcome_model,
      treat = "gdm",
      mediator = "TXE",
      boot = TRUE,
      sims = 100
    )
    
    data.frame(
      transcript_id = txids[tx],
      estimate      = mediate_result$d0,
      lower_95CI    = mediate_result$d0.ci[1],
      upper_95CI    = mediate_result$d0.ci[2],
      pvalue        = mediate_result$d0.p
    )
    
  }, error = function(e) {
    
    warning(paste("Error in transcript", txids[tx], ":", e$message))
    
    data.frame(
      transcript_id = txids[tx],
      estimate      = NA,
      lower_95CI    = NA,
      upper_95CI    = NA,
      pvalue        = NA
    )
  })
}

row.names(Gen3G.mediation.gdm.sig.a.missing) <- NULL
Gen3G.mediation.gdm.sig.a.missing$estimate_abs <- 
  abs(Gen3G.mediation.gdm.sig.a.missing$estimate)

Gen3G.mediation.gdm.sig.a.missing <- 
  left_join(Gen3G.mediation.gdm.sig.a.missing, tx2g.assembly)

Gen3G.mediation.gdm.sig.a.missing$annotation <- "lr-assembly"

Gen3G.mediation.gdm.sig.a.union <- rbind(Gen3G.mediation.gdm.sig.a[,keep_cols],
                                         Gen3G.mediation.gdm.sig.a.union[,keep_cols],
                                         Gen3G.mediation.gdm.sig.a.missing[,keep_cols])
Gen3G.mediation.gdm.sig.a.union <- Gen3G.mediation.gdm.sig.a.union[Gen3G.mediation.gdm.sig.a.union$transcript_id%in%keep_ids,]



# Left join with annotated data to get gene symbols and functions
GUSTO.mediation.gdm.sig.a.union <- GUSTO.mediation.gdm.sig.a.union %>% left_join(union_mapping)
Gen3G.mediation.gdm.sig.a.union <- Gen3G.mediation.gdm.sig.a.union %>% left_join(union_mapping)

GUSTO.mediation.gdm.sig.a.union$label <- factor(
  GUSTO.mediation.gdm.sig.a.union$label,
  levels = unique(GUSTO.mediation.gdm.sig.a.union$label)
)

Gen3G.mediation.gdm.sig.a.union$label <- factor(
  Gen3G.mediation.gdm.sig.a.union$label,
  levels = unique(Gen3G.mediation.gdm.sig.a.union$label)
)

GUSTO.mediation.gdm.sig.a.union$Cohort <- "GUSTO"
Gen3G.mediation.gdm.sig.a.union$Cohort <- "Gen3G"
mediation.gdm.sig.a.union <- rbind(GUSTO.mediation.gdm.sig.a.union,
                                   Gen3G.mediation.gdm.sig.a.union)

# Copy for plotting
mediation.gdm.sig.a.union.plot <- mediation.gdm.sig.a.union

mediation.gdm.sig.a.union.plot$label <- factor(
  mediation.gdm.sig.a.union.plot$label,
  levels = unique(mediation.gdm.sig.a.union.plot$label)
)

mediation.gdm.sig.a.union.plot$Cohort <- factor(
  mediation.gdm.sig.a.union.plot$Cohort,
  levels = c("GUSTO", "Gen3G")
)

# Add significance indicator
mediation.gdm.sig.a.union.plot$significant <- (mediation.gdm.sig.a.union.plot$pvalue <= 0.05)

# Convert to shape numbers
mediation.gdm.sig.a.union.plot$point_shape <- ifelse(mediation.gdm.sig.a.union.plot$significant, 19, 21)

# Create color mapping
function_colors <- c("Cell Organization" = "#F4A6C3",  # pink
                     "Metabolism" = "#F5A623",          # orange
                     "Development" = "#50B8E8",         # light blue
                     "Immune" = "#C5A888",              # tan/beige
                     "Regulation" = "#90D494",          # light green
                     "Unknown" = "#808080")             # grey

# Set fill color - use Function color for significant, white for non-significant
mediation.gdm.sig.a.union.plot$point_fill <- ifelse(mediation.gdm.sig.a.union.plot$significant, 
                                                    function_colors[as.character(mediation.gdm.sig.a.union.plot$Function)], 
                                                    "white")

library(ggtext)

mediation.gdm.sig.a.union.plot$label_formatted <- sapply(as.character(mediation.gdm.sig.a.union.plot$label), function(x) {
  parts <- strsplit(x, "\n")[[1]]
  gene <- parts[1]
  transcript <- parts[2]
  paste0("**", gene, "**<br>", transcript)
})

################################################################################
# Figure 4A: mediation of significant transcripts, union of cohorts
################################################################################

fig4a <- ggplot(mediation.gdm.sig.a.union.plot, aes(x = estimate, y = label_formatted, color = Function)) + 
  geom_errorbarh(aes(xmax = upper_95CI, xmin = lower_95CI), size = 1, height = .2) +
  geom_point(aes(shape = point_shape, fill = point_fill), size = 3) +
  scale_shape_identity() +
  scale_fill_identity() +
  facet_wrap(~Cohort) +
  theme_bw() +
  scale_color_manual(values = function_colors) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white", color = "black"),
        axis.text.y = element_markdown()) +
  ylab("") +
  xlab("Total indirect effect") +
  ggtitle("GDM-birth weight mediation via placental isoform expression") +
  labs(subtitle = "Union of top 10 significant mediators across cohorts\nFilled points are significant at p < 0.05")
fig4a

ggsave("Figures/fig4a.pdf", fig4a, width = 8.39, height = 6.45, units = "in", dpi = 300)


################################################################################
# CSH1 Transcript Structure Visualization (Figure 4B)
################################################################################

# -----------------------------
# 1. Import protein domains
# -----------------------------
import_domtblout_as_granges <- function(domtbl_file, sqanti_table, evalue_thresh = 1e-5, merge_same_domains = TRUE) {
  # Read table
  dom_df <- read.table(domtbl_file, header = TRUE, stringsAsFactors = FALSE)
  
  # Filter high-confidence domains
  dom_df <- dplyr::filter(dom_df, domain_i.Evalue <= evalue_thresh)
  
  # Keep only necessary columns
  dom_df <- dplyr::select(dom_df, target_name, query_name, ali_from, ali_to)
  
  # Join with SQANTI3 CDS info
  dom_df <- dplyr::left_join(
    dom_df,
    dplyr::select(sqanti_table, isoform, CDS_start, CDS_genomic_start, strand, chrom),
    by = c("target_name" = "isoform")
  )
  
  # Filter rows with missing CDS
  dom_df <- dplyr::filter(dom_df, !is.na(CDS_start))
  
  # Map protein coordinates to genomic coordinates
  dom_df <- dplyr::rowwise(dom_df)
  dom_df <- dplyr::mutate(dom_df,
                          domain_tx_start = CDS_start + (as.numeric(ali_from) - 1) * 3,
                          domain_tx_end   = CDS_start + (as.numeric(ali_to)) * 3 - 1,
                          domain_genomic_start = ifelse(strand == "+",
                                                        CDS_genomic_start + (domain_tx_start - CDS_start),
                                                        CDS_genomic_start - (domain_tx_end - CDS_start)),
                          domain_genomic_end   = ifelse(strand == "+",
                                                        CDS_genomic_start + (domain_tx_end - CDS_start),
                                                        CDS_genomic_start - (domain_tx_start))
  )
  dom_df <- dplyr::ungroup(dom_df)
  
  # Ensure start <= end for GRanges
  dom_df <- dom_df %>%
    dplyr::mutate(
      g_start = pmin(domain_genomic_start, domain_genomic_end),
      g_end   = pmax(domain_genomic_start, domain_genomic_end)
    )
  
  # Convert to GRanges - keep domain information
  gr <- GenomicRanges::GRanges(
    seqnames = dom_df$chrom,
    ranges   = IRanges::IRanges(start = dom_df$g_start,
                                end   = dom_df$g_end),
    strand   = dom_df$strand,
    transcript = dom_df$target_name,
    domain     = dom_df$query_name
  )
  
  # Merge overlapping domains if they have the same name
  if (merge_same_domains) {
    # Get unique domain types
    unique_domains <- unique(S4Vectors::mcols(gr)$domain)
    
    # For each domain type, merge overlapping ranges
    gr_list <- lapply(unique_domains, function(dom) {
      gr_subset <- gr[S4Vectors::mcols(gr)$domain == dom]
      gr_merged <- GenomicRanges::reduce(gr_subset)
      S4Vectors::mcols(gr_merged)$domains <- dom
      return(gr_merged)
    })
    
    # Combine all merged domain ranges
    gr <- do.call(c, gr_list)
  } else {
    # Keep domain type in metadata without merging
    S4Vectors::mcols(gr)$domains <- S4Vectors::mcols(gr)$domain
  }
  
  return(gr)
}

# -----------------------------
# 2. Map GRanges to dataframe for ggtranscript
# -----------------------------
map_domains_to_genome <- function(domain_gr) {
  # Each unique domain gets its own track
  dom_df <- data.frame(
    xstart  = start(domain_gr),
    xend    = end(domain_gr),
    track   = mcols(domain_gr)$domains,  # Use domain name as track
    feature = mcols(domain_gr)$domains   # Use domain name as feature for coloring
  )
  dom_df
}

# -----------------------------
# 3. Plot transcripts + domain tracks
# -----------------------------
tx.plot <- function(gtf.df, gene.name, title, domain.gr = NULL, highlight_transcripts = NULL,subtitle=NULL) {
  
  # 1. Extract all exons for the gene
  dat_all <- gtf.df[gtf.df$gene_id == gene.name & gtf.df$type == "exon", ]
  if(nrow(dat_all) == 0) stop("No exons found for gene: ", gene.name)
  
  # 2. Compute full gene range BEFORE subsetting
  full_range <- range(dat_all$start, dat_all$end)
  
  # 3. Desired class order
  desired_order <- c("FSM", "ISM", "NIC", "NNC",
                     "Reference", "Genic Genomic", "Antisense",
                     "Fusion", "Intergenic", "Genic Intron")
  present_classes <- intersect(desired_order, unique(dat_all$class))
  dat_all$class <- factor(dat_all$class, levels = present_classes)
  
  # 4. Use all transcripts
  dat <- dat_all
  
  # 5. Add highlight indicator
  if (!is.null(highlight_transcripts)) {
    dat$highlight <- dat$transcript_id %in% highlight_transcripts
  } else {
    dat$highlight <- FALSE
  }
  
  # 6. Sort transcript IDs by class, then by transcript_id
  transcript_levels <- dat %>%
    distinct(transcript_id, class) %>%
    arrange(factor(class, levels = present_classes), transcript_id) %>%
    pull(transcript_id)
  transcript_levels <- rev(transcript_levels)
  dat$transcript_id <- factor(dat$transcript_id, levels = transcript_levels)
  
  # 7. Domain track (optional)
  domain.mapped <- NULL
  y_levels <- transcript_levels
  
  if (!is.null(domain.gr)) {
    gene.chr <- unique(gtf.df$seqnames[gtf.df$gene_id == gene.name])
    if (length(gene.chr) != 1) stop("Could not uniquely determine chromosome for gene: ", gene.name)
    
    # Subset domain GRanges to gene chromosome and gene range
    domain.gr.gene <- domain.gr[as.character(seqnames(domain.gr)) == gene.chr &
                                  start(domain.gr) <= full_range[2] &
                                  end(domain.gr) >= full_range[1]]
    
    if (length(domain.gr.gene) > 0) {
      # Convert to data.frame for plotting - now with individual domain types as tracks
      domain.mapped <- map_domains_to_genome(domain.gr.gene)
      
      # Get unique domain types for separate tracks
      if(!is.null(domain.mapped) && nrow(domain.mapped) > 0) {
        unique_domains <- unique(domain.mapped$track)
        # Add domain tracks to y_levels
        y_levels <- c(y_levels, unique_domains)
        
        domain.mapped$track <- factor(domain.mapped$track, levels = unique_domains)
        domain.mapped$feature <- factor(domain.mapped$feature, levels = unique_domains)
      }
    }
  }
  
  # 8. Create axis text formatting (bold for highlighted transcripts)
  axis_text_face <- ifelse(y_levels %in% highlight_transcripts, "bold", "plain")
  names(axis_text_face) <- y_levels
  
  # 9. Split data for different line widths
  dat_normal <- dat %>% filter(!highlight)
  dat_highlight <- dat %>% filter(highlight)
  
  # 10. Plot
  p <- ggplot()
  
  # Add normal transcripts first
  if(nrow(dat_normal) > 0) {
    p <- p +
      ggtranscript::geom_range(
        data = dat_normal,
        aes(xstart = start, xend = end, y = transcript_id, fill = class)
      )
  }
  
  # Add highlighted transcripts with thicker lines
  if(nrow(dat_highlight) > 0) {
    p <- p +
      ggtranscript::geom_range(
        data = dat_highlight,
        aes(xstart = start, xend = end, y = transcript_id, fill = class),
        size = 1  # Thicker line for highlighted transcripts
      )
  }
  
  # Add introns for all transcripts
  p <- p +
    ggtranscript::geom_intron(
      data = ggtranscript::to_intron(dat, "transcript_id"),
      aes(xstart = start, xend = end, y = transcript_id, strand = strand),
      arrow = grid::arrow(length = unit(0.05, "inches"))
    ) +
    expand_limits(x = full_range) +
    scale_y_discrete(limits = y_levels) +
    theme_minimal() +
    labs(x = gene.chr,
         subtitle=subtitle,
         title=title) +
    theme(
      panel.grid.major.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 8, face = axis_text_face),
      legend.title = element_blank(),
      legend.position = "top",
      plot.title = element_text(hjust = 0),
      plot.title.position = "plot"
    )
  
  if (!is.null(domain.mapped) && nrow(domain.mapped) > 0) {
    p <- p +
      ggtranscript::geom_range(
        data = domain.mapped,
        aes(xstart = xstart, xend = xend, y = track, fill = feature),
        height = 0.5
      )
  }
  
  return(p)
}

# Create CSH1 gene structure plot
names(classif)[1] <- "isoform"
domain_gr <- import_domtblout_as_granges("Supplementary/pfam.domtblout.txt", classif,merge_same_domains = TRUE)

gtf.assembly <- data.frame(rtracklayer::import("Assembly/ESPRESSO_corrected_SC_filtered.gtf.gz"))
gtf.assembly.CSH1ir <- gtf.assembly[gtf.assembly$transcript_id %in% c("ENST00000453363.7",
                                                                      "ENST00000316193.13",
                                                                      "ENST00000329882.8",
                                                                      "ENST00000558284.1",
                                                                      "ENST00000610991.1",
                                                                      "ENST00000558661.1",
                                                                      "ESPRESSO:chr17:12985:70",
                                                                      "ESPRESSO:chr17:12985:62"), ]
gtf.assembly.CSH1ir$gene_id <- sub("\\..*", "", gtf.assembly.CSH1ir$gene_id)

names(classif)[1] <- "transcript_id"
gtf.assembly.CSH1ir <- left_join(gtf.assembly.CSH1ir,classif[,c("transcript_id","structural_category")])
names(gtf.assembly.CSH1ir)[12] <- "class"

fig4b <- tx.plot(gtf.assembly.CSH1ir, "ENSG00000136488", 
                 "Isoforms of CSH1 (ENSG00000136488) mediating GDM-birth weight effects",
                 domain_gr,
                 highlight_transcripts=c("ESPRESSO:chr17:12985:70",
                                         "ESPRESSO:chr17:12985:62"),
                 subtitle="Significant mediating isoforms in bold") +
  scale_fill_manual(values = c("#66c2a5", "#8da0cb", "grey"))

fig4b

ggsave("Figures/fig4b.pdf", fig4b, width = 6.72, height = 5.56, units = "in", dpi = 300)


################################################################################
# End of Analysis
################################################################################