################################################################################
# Placental Long-Read RNA-seq Assembly Quality Control, Filtering, and Characterization
# Author: Sean T. Bresnahan
# Description: Comprehensive quality control and filtering pipeline for placental
#              transcriptome assembly using SQANTI3 output. Applies multi-layered 
#              filtering based on orthogonal evidence from short-read data, CAGE-seq,
#              poly(A) sites, and artifact detection. Includes post-filtering 
#              characterization, tissue enrichment analysis, and protein validation.
#
# Input Datasets:
#   - Assembly/ESPRESSO_corrected_classification.txt
#       SQANTI3 structural classification of assembled transcripts
#   - Assembly/ESPRESSO_corrected_corrected.gtf
#       GTF annotation file for corrected long-read assembly
#   - Assembly/ESPRESSO_corrected_SC_filtered_classification.txt
#       SQANTI3 classification after structural category filtering
#   - Assembly/ESPRESSO_corrected_SC_filtered.gtf
#       GTF annotation for filtered assembly
#   - Supplementary/GUSTO_SJ.all.filtered.gz
#       Short-read splice junction data from STAR aligner (GUSTO cohort n=200)
#   - Assembly/PLACENTA_ASSEMBLY_corrected_filtered.faa
#       Protein sequences (amino acid FASTA) for assembly
#   - Assembly/PLACENTA_ASSEMBLY_corrected_filtered.fasta
#       Nucleotide sequences (DNA FASTA) for assembly
#   - Assembly/PLACENTA_ASSEMBLY_junctions.txt
#       Splice junction annotations for assembly
#   - Assembly/PLACENTA_ASSEMBLY_corrected.gtf
#       Unfiltered GTF annotation for assembly
#   - /rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v45/gencode.v45.annotation.gtf
#       GENCODE v45 reference annotation
#   - Supplementary/pfam.domtblout.txt
#       Pfam protein domain predictions (HMMER output)
#   - Supplementary/CPC2.txt
#       Coding Potential Calculator 2 (CPC2) predictions
#   - Supplementary/longest_orfs.gff3
#       TransDecoder longest ORF predictions
#   - Supplementary/ESPRESSO_corrected_SC_filtered_blastp.outfmt6
#       BLASTP results against placenta MS/MS database (outfmt6 format)
#   - Supplementary/blastp.outfmt6
#       BLASTP results against full UniProt database (outfmt6 format)
#
# Output Files:
#   Data Files:
#     - Assembly/filter_reasons.txt
#         Table documenting filtering reasons for each removed transcript
#     - Assembly/PLACENTA_ASSEMBLY_classification_filtered.txt
#         Filtered classification table after quality control
#     - Assembly/PLACENTA_ASSEMBLY_corrected_filtered.faa
#         Filtered protein sequences (amino acids)
#     - Assembly/PLACENTA_ASSEMBLY_corrected_filtered.fasta
#         Filtered nucleotide sequences
#     - Assembly/PLACENTA_ASSEMBLY_junctions_filtered.txt
#         Filtered splice junction annotations
#     - Assembly/PLACENTA_ASSEMBLY_corrected_filtered.gtf
#         Filtered GTF annotation
#     - Supplementary/tx2g.RData
#         Transcript-to-gene mapping objects (tx2g.assembly, tx2g.gencode)
#
#   Figures - Quality Control Metrics:
#     - Figures/figs2a.png
#         Full-length read distribution by structural category (boxplot)
#     - Figures/figs2b.png
#         Short-read splice junction support by structural category (boxplot)
#     - Figures/figs2c.png
#         Distance from TSS to CAGE/ATAC peak (histogram with filter thresholds)
#     - Figures/figs2d.png
#         Distance from TTS to poly(A) site (histogram with filter thresholds)
#     - Figures/figs2e.png
#         Distance from TTS to poly(A) motif (histogram with filter thresholds)
#     - Figures/figs2f.png
#         Sequencing artifact detection summary (faceted bar plot)
#
#   Figures - Filtered Assembly Characterization:
#     - Figures/figs3a.png
#         Transcript length distribution post-filtering (boxplot)
#     - Figures/figs3b.png
#         Exon count distribution post-filtering (stacked bar chart)
#     - Figures/figs3c.png
#         Known vs novel gene and isoform classification (boxplot)
#     - Figures/figs3d.png
#         Tissue-specific enrichment of top expressed transcripts (jitter plot)
#
#   Figures - Isoform Complexity and Gene Structure:
#     - Figures/fig2a.pdf
#         Isoform complexity comparison (scatter plot with marginal histograms)
#     - Figures/fig2b.pdf
#         CSH1 gene structure visualization with protein domains
#     - Figures/fig2c.pdf
#         CDS support from placenta MS/MS and UniProt (grouped bar plot)
#
#   Figures - Transcript Features:
#     - Figures/figs4a.png
#         Transcript abundance by novelty category (boxplot)
#     - Figures/figs4b.png
#         Transcript length by novelty category (boxplot)
#     - Figures/figs4c.png
#         Exon count by novelty category (boxplot)
#     - Figures/figs4d.png
#         Correlation between transcripts per gene and exons (all transcripts)
#     - Figures/figs4e.png
#         Correlation between transcripts per gene and exons (highly-expressed)
#
#   Figures - Coding Potential:
#     - Figures/figs5a.png
#         Percentage of isoforms with ORFs by structural category (bar plot)
#     - Figures/figs5b.png
#         ORF length distribution by structural category (boxplot)
#     - Figures/figs5c.png
#         Coding probability distribution by structural category (boxplot)
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
library(seqinr)        # Sequence analysis
library(TissueEnrich)  # Tissue-specific gene enrichment analysis
library(GSEABase)      # Gene set enrichment analysis base functions
library(data.table)    # Enhanced data manipulation
library(enrichR)       # Gene enrichment analysis interface
library(org.Hs.eg.db)  # Human genome annotation database
library(ggtranscript)  # Transcript structure visualization
library(Biostrings)
library(parallel)

################################################################################
# Global Configuration
################################################################################

# Define color palette for structural categories (consistent with previous analysis)
myPalette <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#ffd92f",
               "#e5c494", "#c87774", "#d3adaf", "#b3b3b3")

################################################################################
# Data Import and Preprocessing
################################################################################

# Import SQANTI3 classification results
classif <- read.table("Assembly/ESPRESSO_corrected_classification.txt", header = TRUE)

# Define structural category hierarchy and factor levels
# Establishes order from most to least reliable transcript categories
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

# Rename categories for cleaner visualization and consistent terminology
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

# Handle missing values in minimum coverage data
# Set NA values to 0 for downstream filtering operations
classif[is.na(classif$min_cov), "min_cov"] <- 0

################################################################################
# Short-Read Splice Junction Support Analysis
################################################################################

# Import assembly GTF for transcript structure analysis
GTF <- import("Assembly/ESPRESSO_corrected_corrected.gtf")

# Import short-read splice junction data from STAR aligner
# This provides orthogonal validation of splice sites
SJ.tab <- read.table(gzfile("Supplementary/GUSTO_SJ.all.filtered.gz"), header = FALSE, sep = "\t")

# Assign column names based on STAR SJ.out.tab format
names(SJ.tab) <- c("chr", "start", "end", "strand", "motif", "annotation", 
                   "uniq_reads", "multi_reads", "max_overhang")

# Convert STAR strand encoding to standard notation
# STAR uses: 0=unstranded, 1=+, 2=-
SJ.tab[SJ.tab$strand == 0, "strand"] <- "."
SJ.tab[SJ.tab$strand == 1, "strand"] <- "+"
SJ.tab[SJ.tab$strand == 2, "strand"] <- "-"

# Calculate total read support per splice junction
SJ.tab$total_reads <- SJ.tab$uniq_reads + SJ.tab$multi_reads

# Convert to GenomicRanges object for efficient overlap operations
SJ.tab <- makeGRangesFromDataFrame(SJ.tab, keep.extra.columns = TRUE)

# Aggregate duplicate splice junctions (can occur from multiple samples)
df <- as.data.frame(SJ.tab)
df_summary <- df %>%
  group_by(seqnames, start, end, strand, motif, annotation) %>%
  dplyr::summarise(
    uniq_reads = sum(uniq_reads),
    multi_reads = sum(multi_reads),
    max_overhang = max(max_overhang),
    .groups = "drop"
  )

# Recreate GenomicRanges object with aggregated data
SJ.tab <- makeGRangesFromDataFrame(
  df_summary,
  keep.extra.columns = TRUE,
  seqnames.field = "seqnames",
  start.field = "start",
  end.field = "end",
  strand.field = "strand"
)

# Recalculate total reads after aggregation
SJ.tab$total_reads <- SJ.tab$uniq_reads + SJ.tab$multi_reads

################################################################################
# TSS Ratio Calculation from Short-Read Data
################################################################################

# Extract exon annotations from assembly GTF
gtf_exons <- GTF[GTF$type == "exon"]
gtf_exons <- data.frame(gtf_exons)

# Number exons within each transcript for ratio calculation
gtf_exons <- gtf_exons %>%
  dplyr::group_by(transcript_id) %>%
  dplyr::arrange(seqnames, start, end, .by_group = TRUE) %>%
  dplyr::mutate(exon_number = row_number()) %>%
  ungroup()

# Convert back to GenomicRanges for overlap analysis
gtf_exons <- makeGRangesFromDataFrame(gtf_exons, keep.extra.columns = TRUE)

# Find overlaps between exons and splice junctions
# This identifies which splice junctions support each exon
hits <- findOverlaps(gtf_exons, SJ.tab)
exon.hits <- gtf_exons[queryHits(hits)]
SJ.hits <- SJ.tab[subjectHits(hits)]

# Transfer read counts from splice junctions to exons
exon.hits$total_reads <- SJ.hits$total_reads

# Summarize read counts per exon per transcript
summarized_exons <- as.data.frame(exon.hits) %>%
  dplyr::group_by(transcript_id, exon_number) %>%
  dplyr::summarise(total_reads_sum = sum(total_reads), .groups = "drop")

# Calculate TSS ratio: first exon reads / last exon reads
# Lower ratios may indicate degradation or incomplete transcripts
tss_ratio <- summarized_exons %>%
  dplyr::group_by(transcript_id) %>%
  dplyr::summarise(
    first_exon_reads = total_reads_sum[exon_number == min(exon_number)],
    last_exon_reads = total_reads_sum[exon_number == max(exon_number)],
    TSS_ratio = first_exon_reads / last_exon_reads
  )

# Standardize column name for merging
names(tss_ratio)[1] <- "isoform"

# Add TSS ratio to classification data
classif <- left_join(classif, tss_ratio[, c(1, 4)])

################################################################################
# Quality Control Filter Functions
################################################################################

# Generic filtering function that applies a condition and tracks filtered transcripts
filter.classif <- function(data, func, label) {
  # Apply filter condition to retain high-quality transcripts
  data.filt <- data %>% filter(eval(parse(text = func)))
  
  # Identify transcripts that failed the filter
  out <- data.frame(isoform = setdiff(data$isoform, data.filt$isoform))
  out <- left_join(out, data[, c("isoform", "structural_category")], by = "isoform")
  out$reason <- label
  
  return(out)
}

# Logical OR operation for combining filter results
# Transcripts failing either filter are flagged
filter.or <- function(a, b) {
  int <- intersect(a$isoform, b$isoform)
  return(rbind(a[a$isoform %in% int, ], b[b$isoform %in% int, ]))
}

# Logical AND operation for combining filter results  
# Transcripts failing any filter are flagged
filter.and <- function(a, b) {
  int <- union(a$isoform, b$isoform)
  return(rbind(a[a$isoform %in% int, ], b[b$isoform %in% int, ]))
}

################################################################################
# Orthogonal Evidence-Based Filtering
################################################################################

# Initialize data frame to track filter results
dat.figs2f <- classif[, c("isoform", "structural_category")]

## 3' End Support (TTS - Transcription Termination Site)
# Multiple lines of evidence for proper 3' end annotation

# Filter 1: Distance to poly(A) site annotation
# Transcripts should terminate near known poly(A) sites (±100bp)
filt.dist_to_polyA_site <- filter.classif(classif, 
                                          "dist_to_polyA_site > -100 & dist_to_polyA_site < 100", 
                                          "dist_to_polyA_site")

# Filter 2: Distance to poly(A) motif  
# Transcripts should have poly(A) motifs upstream of cleavage site (-41 to 0bp)
filt.polyA_dist <- filter.classif(classif, 
                                  "polyA_dist > -41 & polyA_dist < 1", 
                                  "polyA_dist")

# Filter 3: Difference to annotated TTS
# Transcripts should end near reference transcript termination sites (±100bp)
filt.diff_to_TTS <- filter.classif(classif, 
                                   "diff_to_TTS > -100 & diff_to_TTS < 100", 
                                   "diff_to_TTS")

# Combine 3' end filters (transcript passes if any condition is met)
filt.3prime <- filter.or(filter.or(filt.dist_to_polyA_site, filt.polyA_dist), filt.diff_to_TTS)

# Mark transcripts with 3' end support
dat.figs2f$TTS <- FALSE
dat.figs2f[dat.figs2f$isoform %in% filt.3prime$isoform, "TTS"] <- TRUE

## 5' End Support (TSS - Transcription Start Site)
# Evidence for proper transcription initiation

# Filter 1: Distance to CAGE peak
# CAGE-seq identifies authentic transcription start sites (±5kb tolerance)
filt.dist_to_CAGE_peak <- filter.classif(classif, 
                                         "dist_to_CAGE_peak > -5000 & dist_to_CAGE_peak < 5000", 
                                         "dist_to_CAGE_peak")

# Filter 2: Difference to annotated TSS
# Transcripts should start near reference transcript start sites (±100bp)
filt.diff_to_TSS <- filter.classif(classif, 
                                   "diff_to_TSS > -100 & diff_to_TSS < 100", 
                                   "diff_to_TSS")

# Combine 5' end filters
filt.5prime <- filter.or(filt.dist_to_CAGE_peak, filt.diff_to_TSS)

# Mark transcripts with 5' end support
dat.figs2f$TSS <- FALSE
dat.figs2f[dat.figs2f$isoform %in% filt.5prime$isoform, "TSS"] <- TRUE

## Splice Junction Support
# Short-read validation of splice sites (minimum 3 reads per junction)
filt.min_cov <- filter.classif(classif, "min_cov > 2", "min_cov")

dat.figs2f$SJ <- FALSE
dat.figs2f[dat.figs2f$isoform %in% filt.min_cov$isoform, "SJ"] <- TRUE

## Full-Length Read Support
# Long-read evidence for complete transcript structure (minimum 3 reads)
filt.FL <- filter.classif(classif, "FL > 2", "FL")

dat.figs2f$FL <- FALSE
dat.figs2f[dat.figs2f$isoform %in% filt.FL$isoform, "FL"] <- TRUE

################################################################################
# Artifact Detection and Filtering
################################################################################

## Poly(A) Internal Priming Artifact Detection
# Internal priming occurs when oligo(dT) primes within A-rich regions
# Exclude FSM transcripts (reference-validated) from this filter
filt.perc_A_downstream_TTS <- filter.classif(classif[!classif$structural_category == "FSM", ],
                                             "perc_A_downstream_TTS > -1 & perc_A_downstream_TTS < 80", 
                                             "perc_A_downstream_TTS")

# Mark transcripts with poly(A) intrapriming artifacts
dat.figs2f$perc_A_downstream_TTS <- FALSE
dat.figs2f[dat.figs2f$isoform %in% filt.perc_A_downstream_TTS$isoform, "perc_A_downstream_TTS"] <- TRUE

## RT-Switching Artifact Detection
# Template switching during reverse transcription creates false splice junctions
# Apply only to novel transcript categories (FSM, ISM, NIC are reference-validated)
filt.RTS_stage <- filter.classif(classif[!classif$structural_category %in% c("FSM", "ISM", "NIC"), ],
                                 "RTS_stage == FALSE", 
                                 "RTS_stage")

dat.figs2f$RTS_stage <- FALSE
dat.figs2f[dat.figs2f$isoform %in% filt.RTS_stage$isoform, "RTS_stage"] <- TRUE

## Degradation Product Detection
# Identify transcripts likely to be degradation products or NMD targets

# Filter 1: Predicted nonsense-mediated decay (NMD)
# Transcripts with premature stop codons are NMD targets
filt.NMD <- filter.classif(classif[!classif$structural_category %in% c("FSM", "ISM", "NIC"), ],
                           "predicted_NMD == FALSE", 
                           "predicted_NMD")

# Filter 2: TSS ratio indicates degradation
# Low TSS ratios suggest 5' degradation
filt.TSS_ratio <- filter.classif(classif[!classif$structural_category %in% c("FSM", "ISM", "NIC"), ],
                                 "TSS_ratio > 0.5", 
                                 "TSS_ratio")

# Combine degradation filters
filt.predicted_NMD <- filter.or(filt.NMD, filt.TSS_ratio)

# Mark transcripts as potential degradation products
dat.figs2f$Decay <- FALSE
dat.figs2f[dat.figs2f$isoform %in% filt.predicted_NMD$isoform, "Decay"] <- TRUE

################################################################################
# Summary Statistics for Figure S2F
################################################################################

# Transform data from wide to long format for aggregation
# This prepares the data for calculating percentages by structural category
dat.figs2f.long <- dat.figs2f[, -1] %>%
  pivot_longer(
    cols = c(TTS, TSS, SJ, FL, perc_A_downstream_TTS, RTS_stage, Decay),
    names_to = "type",
    values_to = "value"
  )

# Convert logical values to integers for aggregation
dat.figs2f.long$value <- as.integer(dat.figs2f.long$value)

# Calculate mean proportion for each combination of structural category and filter type
# This gives the percentage of transcripts passing each filter within each category
dat.figs2f.summary <- aggregate(value ~ structural_category + type, 
                                data = dat.figs2f.long, 
                                FUN = mean)

# Rename filter types for cleaner visualization labels
dat.figs2f.summary$type <- revalue(dat.figs2f.summary$type, 
                                   c("FL" = "Full-length reads",
                                     "SJ" = "Splice-junction reads", 
                                     "perc_A_downstream_TTS" = "Poly(A) intrapriming",
                                     "RTS_stage" = "RT switching",
                                     "Decay" = "NMD or degradation product"))

################################################################################
# Generate Filtered Assembly Files
################################################################################

# Note: Some filter variables appear to be missing from the provided code
# This section would typically combine all filters to create the final filtered set

## Compile all filtering reasons
# This tracks why each transcript was filtered for quality control
filt.cov <- filter.and(filt.FL, filt.min_cov)
reasons <- rbind(
  filt.cov,
  filt.5prime,
  filt.3prime,
  filt.perc_A_downstream_TTS,
  filt.RTS_stage,
  filt.predicted_NMD
)

# Write filtering report
# write.table(reasons, "Assembly/filter_reasons.txt",
#             sep = "\t", quote = FALSE, row.names = FALSE)

# Create filtered classification table
classif.filtered <- classif[!classif$isoform %in% reasons$isoform, ]
# write.table(classif.filtered, "Assembly/PLACENTA_ASSEMBLY_classification_filtered.txt",
#             sep = "\t", quote = FALSE, row.names = FALSE)

################################################################################
# Filter Supporting Files
################################################################################

## Filter protein sequences (amino acid FASTA)
faa <- read.fasta(file = "Assembly/PLACENTA_ASSEMBLY_corrected_filtered.faa", 
                  seqtype = "AA", as.string = TRUE, set.attributes = FALSE)

# Extract transcript IDs from FASTA headers and filter
faa <- faa[sapply(strsplit(names(faa), "\t"), `[[`, 1) %in% classif.filtered$isoform]
write.fasta(sequences = faa, names = names(faa), 
            file.out = paste0(DIR_OUT, "/", OUT_PREFIX, ".faa"))
rm(faa)

## Filter nucleotide sequences (DNA FASTA)
fasta <- read.fasta(file = "Assembly/PLACENTA_ASSEMBLY_corrected_filtered.fasta", 
                    seqtype = "DNA", as.string = TRUE, set.attributes = FALSE)
fasta <- fasta[names(fasta) %in% classif.filtered$isoform]
write.fasta(sequences = fasta, names = names(fasta), 
            file.out = paste0(DIR_OUT, "/", OUT_PREFIX, ".fasta"))
rm(fasta)

## Filter splice junction annotations
junc <- read.table("Assembly/PLACENTA_ASSEMBLY_junctions.txt", header = TRUE)
junc <- junc[junc$isoform %in% classif.filtered$isoform, ]
write.table(junc, "Assembly/PLACENTA_ASSEMBLY_junctions_filtered.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

## Filter GTF annotation
gtf <- readGFF("Assembly/PLACENTA_ASSEMBLY_corrected.gtf")
gtf <- gtf[gtf$transcript_id %in% classif.filtered$isoform, ]

# Format GTF attributes for standard compliance
gtf$score <- "."
gtf$phase <- "."
gtf$info <- NA

# Create properly formatted attribute strings
gtf[gtf$type == "transcript", "info"] <- paste0("transcript_id \"",
                                                gtf[gtf$type == "transcript", "transcript_id"],
                                                "\"; gene_id \"",
                                                gtf[gtf$type == "transcript", "gene_id"], "\";")

gtf[gtf$type == "exon", "info"] <- paste0("transcript_id \"",
                                          gtf[gtf$type == "exon", "transcript_id"],
                                          "\"; gene_id \"",
                                          gtf[gtf$type == "exon", "gene_id"], "\";")

# Remove redundant columns
gtf$transcript_id <- NULL
gtf$gene_id <- NULL

# Write filtered GTF
write.table(gtf, "Assembly/PLACENTA_ASSEMBLY_corrected_filtered.gtf", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

################################################################################
# Quality Control Visualizations
################################################################################

## Figure S2A: Full-Length Read Distribution
# Shows support level across structural categories
figs2a <- ggplot(classif, aes(x = structural_category, y = FL, fill = structural_category)) +
  geom_boxplot() +
  ylab("Full-length reads in n=72") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  guides(fill = "none") +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_fill_manual(values = myPalette) +
  ggtitle("Full-length reads") +
  # Add filtering threshold line
  geom_hline(yintercept = 3, color = 'black',
             linetype = "dashed", linewidth = 0.5)

figs2a
ggsave("Figures/figs2a.png", figs2a, width = 4, height = 4, units = "in", dpi = 300)

## Figure S2B: Short-Read Splice Junction Support
# Minimum coverage across all splice junctions per transcript
figs2b <- ggplot(classif, aes(x = structural_category, y = min_cov + 1, fill = structural_category)) +
  geom_boxplot() +
  ylab("Short reads in GUSTO n=200") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  guides(fill = "none") +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_fill_manual(values = myPalette) +
  ggtitle("Minimum splice-junction reads") +
  # Add filtering threshold line
  geom_hline(yintercept = 3, color = 'black',
             linetype = "dashed", linewidth = 0.5)

figs2b
ggsave("Figures/figs2b.png", figs2b, width = 4, height = 4, units = "in", dpi = 300)

## Figure S2C: TSS Validation by CAGE/ATAC-seq
# Distance from transcript start to experimentally validated start sites
figs2c <- ggplot(classif, aes(x = dist_to_CAGE_peak, fill = structural_category)) +
  geom_histogram(aes(y = after_stat(count / sum(count))), 
                 binwidth = 100, position = "stack") +
  xlab("bp") + 
  ylab("% of isoforms") +
  theme_minimal() + 
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none") +
  scale_fill_manual(values = myPalette) +
  ggtitle("Distance from TSS to\nCAGE/ATAC peak") +
  xlim(-10000, 10500) +
  # Add filtering threshold lines
  geom_vline(xintercept = -5000, color = 'black',
             linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = 5000, color = 'black',
             linetype = "dashed", linewidth = 0.5)

figs2c
ggsave("Figures/figs2c.png", figs2c, width = 3, height = 3, units = "in", dpi = 300)

## Figure S2D: Poly(A) Site Validation  
# Distance from transcript end to annotated polyadenylation sites
figs2d <- ggplot(classif, aes(x = dist_to_polyA_site, fill = structural_category)) +
  geom_histogram(aes(y = after_stat(count / sum(count))), 
                 binwidth = 1, position = "stack") +
  xlab("bp") + 
  ylab("% of isoforms") +
  theme_minimal() + 
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none") +
  scale_fill_manual(values = myPalette) +
  ggtitle("Distance from TTS to\nPoly(A) site") +
  # Add filtering threshold lines
  geom_vline(xintercept = -100, color = 'black',
             linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = 100, color = 'black',
             linetype = "dashed", linewidth = 0.5)

figs2d
ggsave("Figures/figs2d.png", figs2d, width = 3, height = 3, units = "in", dpi = 300)

## Figure S2E: Poly(A) Motif Distance Analysis
# Distance from transcript end to poly(A) signal motifs
figs2e <- ggplot(classif, aes(x = polyA_dist, fill = structural_category)) +
  geom_histogram(aes(y = after_stat(count / sum(count))), 
                 binwidth = 1, position = "stack") +
  xlab("bp") + 
  ylab("% of isoforms") +
  theme_minimal() + 
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none") +
  scale_fill_manual(values = myPalette) +
  ggtitle("Distance from TTS to\nPoly(A) motif") +
  # Add filtering threshold lines
  geom_vline(xintercept = -41, color = 'black',
             linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = 'black',
             linetype = "dashed", linewidth = 0.5)

figs2e
ggsave("Figures/figs2e.png", figs2e, width = 3, height = 3, units = "in", dpi = 300)

## Figure S2F: Sequencing Artifact Detection Summary
# Shows the percentage of transcripts affected by each type of artifact by structural category
figs2f <- ggplot(dat.figs2f.summary[dat.figs2f.summary$type %in% c("Poly(A) intrapriming", 
                                                                   "RT switching",
                                                                   "NMD or degradation product"), ], 
                 aes(x = structural_category, y = value, fill = structural_category)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ type) +
  ylab("% of isoforms") +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  guides(fill = "none") +
  scale_fill_manual(values = myPalette) +
  ggtitle("Sequencing artifacts")

figs2f
ggsave("Figures/figs2f.png", figs2f, width = 6, height = 3, units = "in", dpi = 300)

################################################################################
# Filtered Assembly Characterization and Visualization
################################################################################

## Figure S3A: Transcript Length Distribution (Post-Filtering)
# Shows length characteristics of high-quality transcripts after QC filtering
figs3a <- ggplot(classif.filtered, aes(y = length, x = structural_category, fill = structural_category)) +
  geom_boxplot() + 
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_fill_manual(values = myPalette) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "none") +
  labs(title = "Isoform length",
       y = "Length (bp)") 

figs3a
ggsave("Figures/figs3a.png", figs3a, width = 3, height = 3, units = "in", dpi = 300)

## Figure S3B: Exon Count Distribution (Post-Filtering)
# Distribution of exon counts across structural categories after filtering
figs3b <- ggplot(classif.filtered, aes(x = exons, fill = structural_category)) +
  geom_bar(stat = "count", position = "stack") +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values = myPalette) +
  theme_minimal() + 
  theme(legend.position = "none") +
  labs(title = "Exons per isoform by structural category",
       x = "Exons", 
       y = "Isoforms") +
  xlim(0, 25)  # Focus on reasonable exon count range

figs3b
ggsave("Figures/figs3b.png", figs3b, width = 4, height = 4, units = "in", dpi = 300)

## Figure S3C: Known vs Novel Gene and Isoform Classification
# Categorizes transcripts by novelty at both gene and isoform levels
figs3c.dat <- classif.filtered[, c("associated_gene", "structural_category")]

# Classify genes as known (ENSEMBL) or novel (coordinate-based)
figs3c.dat$gene_class <- "Known"
figs3c.dat[grep("novel", figs3c.dat$associated_gene), "gene_class"] <- "Novel"

# Classify isoforms as known (FSM = Full Splice Match) or novel
figs3c.dat$isoform_class <- "Novel"
figs3c.dat[figs3c.dat$structural_category == "FSM", "isoform_class"] <- "Known"

# Remove structural category for cleaner analysis
figs3c.dat$structural_category <- NULL

# Count combinations and filter invalid combinations
figs3c.dat <- data.frame(table(figs3c.dat))
figs3c.dat$keep <- TRUE

# Remove impossible combination: novel gene with known isoform
figs3c.dat[figs3c.dat$gene_class == "Novel" & figs3c.dat$isoform_class == "Known", "keep"] <- FALSE

# Remove empty combinations
figs3c.dat[figs3c.dat$gene_class == "Known" & figs3c.dat$isoform_class == "Known" & figs3c.dat$Freq == 0, "keep"] <- FALSE

figs3c.dat <- figs3c.dat[figs3c.dat$keep == TRUE, ]

# Note: The original code references 'g7.dat' which appears to be a typo for 'figs3c.dat'
figs3c <- ggplot(figs3c.dat, aes(y = log(Freq + 1), x = gene_class, fill = isoform_class)) +
  geom_boxplot() + 
  scale_fill_manual(values = c("grey", "white")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "top") +
  labs(title = "Known vs novel",
       y = "Isoforms (x10³)", 
       x = "Gene",
       fill = "Isoform")

figs3c
ggsave("Figures/figs3c.png", figs3c, width = 3, height = 3, units = "in", dpi = 300)

################################################################################
# Figure S3D: Tissue-Specific Enrichment Analysis
################################################################################

# Regenerate transcript-to-gene mapping for filtered assembly
tx2g.assembly <- data.frame(import("Assembly/ESPRESSO_corrected_SC_filtered.gtf"))
tx2g.assembly <- tx2g.assembly[!duplicated(tx2g.assembly), ]
tx2g.assembly$score <- NULL
tx2g.assembly$phase <- NULL

# Handle novel transcripts without gene assignments
novelTx.assembly <- unique(tx2g.assembly[tx2g.assembly$gene_id == "NA", "transcript_id"])

# Use SQANTI3 associated_gene assignments for novel transcripts when available
novelTx.ascg <- classif[classif$isoform %in% novelTx.assembly, 
                        c("isoform", "associated_gene")]
novelTx.ascg <- novelTx.ascg[grep("ENSG", novelTx.ascg$associated_gene), ]
names(novelTx.ascg)[1] <- "transcript_id"

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

# Make tx2g for gencode
tx2g.gencode <- data.frame(import("/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v45/gencode.v45.annotation.gtf"))
tx2g.gencode <- tx2g.gencode[, c("transcript_id", "gene_id")]
tx2g.gencode <- tx2g.gencode[!duplicated(tx2g.gencode), ]

# Save tx2g objects for later use
save(list=c("tx2g.assembly","tx2g.gencode"),file="Supplementary/tx2g.RData")

# Function to convert ENSEMBL IDs to gene symbols
ensembl_to_symbol <- function(ensembl_ids) {
  gene_symbols <- mapIds(
    org.Hs.eg.db,
    keys = ensembl_ids,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  # Remove NAs and create mapping
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

# Prepare gene universe for enrichment analysis
tx2g.assembly$gene_id <- sapply(strsplit(tx2g.assembly$gene_id, '[.]'),
                                function(x) x[1])
universe_assembly <- unique(tx2g.assembly$gene_id)
universe_assembly_result <- ensembl_to_symbol(universe_assembly)
symbols_universe_assembly <- universe_assembly_result$symbols
mapping_universe_assembly <- universe_assembly_result$mapping

# Rank genes by expression (full-length read count) for enrichment analysis
gexp <- classif.filtered[rev(order(classif.filtered$FL)), "associated_gene"]

# Initialize results data frame for tissue enrichment analysis
results.df <- data.frame(matrix(ncol = 8, nrow = 0))
names(results.df) <- c("Log10PValue", "Tissue.Specific.Genes", "fold.change", "samples",
                       "Tissue", "OddsRatio", "p", "Ntop")

# Define gene set sizes to test (top 100 to 4000 genes)
sets <- seq.int(100, 4000, 100)

# Perform tissue enrichment analysis across different top gene set sizes
for(set in sets) {
  # Extract top N expressed genes
  Xset <- gexp[1:set]
  Xset <- sapply(strsplit(Xset, '[.]'), function(x) x[1])
  symbols_Xset <- ensembl_to_symbol(Xset)[["mapping"]][["gene_symbol"]]
  
  # Create gene sets for TissueEnrich analysis
  gs <- GeneSet(geneIds = unique(symbols_Xset), 
                organism = 'Homo Sapiens',
                geneIdType = SymbolIdentifier())
  
  bs <- GeneSet(geneIds = unique(symbols_universe_assembly), 
                organism = 'Homo Sapiens',
                geneIdType = SymbolIdentifier())
  
  # Perform tissue enrichment
  results <- teEnrichment(gs, backgroundGenes = bs)
  enrich_table <- data.frame(results[[1]]@assays@data@listData[[1]])
  
  # Load tissue-specific gene data for odds ratio calculation
  bg_genes <- unique(mapping_universe_assembly$ensembl_id)
  tissue_data <- read.csv("/Users/stbresnahan/Desktop/Manuscripts/Placenta/fig5_SC/EnrichR_Results_MSSupp/tsgs.csv")
  tissue_data$Tissue <- gsub("\\.", " ", tissue_data$Tissue)
  
  # Calculate odds ratios for each tissue
  enrich_table$Tissue <- row.names(enrich_table)
  enrich_table$OddsRatio <- NA
  my_genes <- unique(Xset)
  
  for (i in 1:nrow(enrich_table)) {
    tissue_name <- enrich_table$Tissue[i]
    tissue_genes <- tissue_data[tissue_data$Tissue == tissue_name, ]
    
    # Construct 2x2 contingency table for Fisher's exact test
    a <- sum(my_genes %in% tissue_genes$Gene)                                      # in list & tissue
    b <- sum(my_genes %in% bg_genes & !(my_genes %in% tissue_genes$Gene))         # in list & not tissue
    c <- sum(!(bg_genes %in% my_genes) & (bg_genes %in% tissue_genes$Gene))       # not in list & in tissue
    d <- sum(!(bg_genes %in% my_genes) & !(bg_genes %in% tissue_genes$Gene))      # not in list & not in tissue
    
    contingency_table <- matrix(c(a, b, c, d), nrow = 2)
    ft <- fisher.test(contingency_table)
    enrich_table$OddsRatio[i] <- ft$estimate
  }
  
  # Add metadata and combine results
  enrich_table$p <- 10^-enrich_table$Log10PValue
  enrich_table$Ntop <- set
  results.df <- rbind(results.df, enrich_table)
}

# Adjust p-values and prepare for visualization
results.df$neglog10p <- -log10(results.df$p)
results.df$padj <- p.adjust(results.df$p, method = "fdr")
results.df$neglog10padj <- -log10(results.df$padj)
results.df$Tissue <- as.factor(results.df$Tissue)

# Set up color scheme with placenta highlighted
default_colors <- setNames(rainbow(length(levels(results.df$Tissue))), 
                           levels(results.df$Tissue))
default_colors["Placenta"] <- "black"

# Create tissue enrichment plot
figs3d <- ggplot(results.df[results.df$Ntop %in% c(100:2000), ], 
                 aes(x = Ntop, y = OddsRatio)) +
  geom_boxplot(aes(group = Ntop), color = "grey40", outliers = FALSE) +
  geom_jitter(aes(color = Tissue), width = 15) +
  scale_color_manual(values = default_colors) +
  theme_minimal() +
  xlab("Top expressed transcripts") +
  ylab("Odds ratio") +
  ggtitle("Tissue-specific enrichment of top expressed transcripts") +
  theme(legend.position = "bottom")

figs3d
ggsave("Figures/figs3d.png", figs3d, width = 10.55, height = 5.92, units = "in", dpi = 300)

################################################################################
# Figure 2A: Isoform Complexity Comparison (Filtered Assembly)
################################################################################

# Prepare filtered assembly transcript-to-gene mapping
tx2g.assembly.filtered <- classif.filtered[, c("isoform", "associated_gene")]
names(tx2g.assembly.filtered)[2] <- "gene_id"

# Calculate isoforms per gene for filtered assembly
iso.tab.assembly.filtered <- data.frame(table(tx2g.assembly.filtered$gene_id))
names(iso.tab.assembly.filtered) <- c("gene_id", "Freq.assembly")
iso.tab.assembly.filtered$gene_id <- as.character(iso.tab.assembly.filtered$gene_id)

# Calculate isoforms per gene for gencode
iso.tab.gencode <- data.frame(table(tx2g.gencode$gene_id))
names(iso.tab.gencode) <- c("gene_id", "Freq.gencode")
iso.tab.gencode$gene_id <- as.character(iso.tab.gencode$gene_id)

# Create comprehensive comparison dataset (filtered vs reference)
genes.all <- union(tx2g.assembly.filtered$gene_id, tx2g.gencode$gene_id)
compgene.filtered <- data.frame(gene_id = genes.all)
compgene.filtered <- left_join(compgene.filtered, iso.tab.assembly.filtered, by = "gene_id")
compgene.filtered <- left_join(compgene.filtered, iso.tab.gencode, by = "gene_id")
compgene.filtered[is.na(compgene.filtered)] <- 0

# Classify genes by isoform complexity change
compgene.filtered$set <- "no change"
compgene.filtered[(compgene.filtered$Freq.gencode - compgene.filtered$Freq.assembly) > 0, "set"] <- "reduced"
compgene.filtered[(compgene.filtered$Freq.gencode - compgene.filtered$Freq.assembly) < 0, "set"] <- "expanded"
compgene.filtered$set <- factor(compgene.filtered$set, levels = c("expanded", "reduced", "no change"))
compgene.filtered <- compgene.filtered[order(compgene.filtered$set), ]

# Highlight specific genes of interest (Growth Hormone cluster)
compgene.filtered$col <- ifelse(compgene.filtered$gene_id %in% c("ENSG00000136488.15", "ENSG00000136487.18"), 
                                "highlight", "normal")
highlighted <- compgene.filtered[compgene.filtered$col == "highlight", ]

# Create scatter plot with highlighted genes
fig2a <- ggplot(compgene.filtered, aes(x = log(Freq.gencode + 1), y = log(Freq.assembly + 1))) +
  geom_point(size = 0.6, aes(color = col)) +
  geom_label(data = highlighted,
             aes(label = c("CSH1", "GH2")),
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
  xlab("GENCODEv45 (log+1)") +
  ylab("lr-assembly (log+1)") +
  ggtitle("Isoforms per gene")

# After creating fig1f with marginal histograms
fig2a_with_margins <- ggMarginal(fig2a, type = "histogram", fill = "grey40", color = "white")

# Save as PDF
ggsave("Figures/fig2a.pdf", 
       plot = fig2a_with_margins,
       width = 4.82, 
       height = 4.82, 
       units = "in",
       dpi = 300,
       device = "pdf")

################################################################################
# Figure 2B: CSH1 Gene Structure Visualization
################################################################################

# Prepare combined GTF data for transcript structure plotting
gtf.assembly <- data.frame(import("Assembly/ESPRESSO_corrected_SC_filtered.gtf"))
gtf.assembly$gene_id <- NULL
names(tx2g.assembly.filtered)[1] <- "transcript_id"
gtf.assembly <- left_join(gtf.assembly, tx2g.assembly.filtered)
lj <- classif[, c("isoform", "structural_category")]
names(lj)[1] <- "transcript_id"
gtf.assembly <- left_join(gtf.assembly, lj)
names(gtf.assembly)[12] <- "class"

# Import reference transcripts for comparison
gtf.gencode <- data.frame(import("/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v45/gencode.v45.annotation.gtf"))
gtf.gencode <- gtf.gencode[, names(gtf.gencode) %in% names(gtf.assembly)]
gtf.gencode$class <- "Reference"
gtf.gencode <- gtf.gencode[!gtf.gencode$transcript_id %in% gtf.assembly$transcript_id, ]
genes_with_transcripts <- unique(gtf.gencode$gene_id[gtf.gencode$type == "transcript"])
gtf.gencode <- gtf.gencode[gtf.gencode$type != "gene" | 
                             gtf.gencode$gene_id %in% genes_with_transcripts, ]

# Combine assembly and reference annotations
gtf.assembly <- rbind(gtf.assembly, gtf.gencode)
gtf.assembly$gene_id <- sapply(strsplit(gtf.assembly$gene_id, '[.]'), function(x) x[1])
rm(gtf.gencode)

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
  desired_order <- c("Reference", "FSM", "ISM", "NIC", "NNC",
                     "Genic Genomic", "Antisense",
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
      color="grey",
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

names(classif.filtered)[1] <- "isoform"
domain_gr <- import_domtblout_as_granges("Supplementary/pfam.domtblout.txt", classif.filtered)

gtf.assembly$seqnames <- as.character(gtf.assembly$seqnames)
gtf.assembly <- gtf.assembly[!is.na(gtf.assembly$seqnames),]
gtf.assembly <- gtf.assembly[!is.na(gtf.assembly$gene_id), ]

fig2b <- tx.plot(gtf.assembly, "ENSG00000136488", "CSH1 (ENSG00000136488)",domain_gr) +
  scale_fill_manual(values = c("#66c2a5", "#8da0cb", "#e78ac3","grey"))

fig2b
ggsave("Figures/fig2b.pdf", fig2b, width = 8.65, height = 11, units = "in", dpi = 300)


################################################################################
# Transcript Feature Analysis
################################################################################

names(classif.filtered)[1] <- "transcript_id"

## Classify transcripts by novelty for feature analysis
novelGtx <- classif.filtered[grep("novel", classif.filtered$associated_gene), "transcript_id"]
fsm <- classif.filtered[classif.filtered$structural_category %in% c("FSM", "ISM"), "transcript_id"]
noveltx <- setdiff(classif.filtered[!classif.filtered$structural_category %in% c("FSM", "ISM"), "transcript_id"],
                   novelGtx)

################################################################################
# Figure S4A: Transcript abundance
################################################################################

fl.dat <- data.frame(FL=c(classif.filtered[classif.filtered$transcript_id%in%fsm,"FL"],
                          classif.filtered[classif.filtered$transcript_id%in%noveltx,"FL"],
                          classif.filtered[classif.filtered$transcript_id%in%c(fsm,noveltx),"FL"],
                          classif.filtered[classif.filtered$transcript_id%in%novelGtx,"FL"]),
                     set=c(rep.int("Known",length(classif.filtered[classif.filtered$transcript_id%in%fsm,"FL"])),
                           rep.int("Novel",length(classif.filtered[classif.filtered$transcript_id%in%noveltx,"FL"])),
                           rep.int("Known",length(classif.filtered[classif.filtered$transcript_id%in%c(fsm,noveltx),"FL"])),
                           rep.int("Novel",length(classif.filtered[classif.filtered$transcript_id%in%novelGtx,"FL"]))),
                     comparison=c(rep.int("Within known genes",length(classif.filtered[classif.filtered$transcript_id%in%fsm,"FL"])+length(classif.filtered[classif.filtered$transcript_id%in%noveltx,"FL"])),
                                  rep.int("Known vs novel genes",length(classif.filtered[classif.filtered$transcript_id%in%c(fsm,noveltx),"FL"])+length(classif.filtered[classif.filtered$transcript_id%in%novelGtx,"FL"]))))

figs4a <- ggplot(fl.dat, aes(x = comparison, y = log10(FL), color = set)) +
  geom_boxplot() + 
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Full-length reads (log10)") +
  scale_color_manual(values = c("black", "grey"),
                     name = "Transcripts")

figs4a
ggsave("figs4a.png", figs4a, width = 4, height = 4, units = "in", dpi = 300)

# Statistical testing
wilcox.test(classif.filtered[classif.filtered$transcript_id %in% fsm, "FL"],
            classif.filtered[classif.filtered$transcript_id %in% noveltx, "FL"],
            alternative = "greater")

wilcox.test(classif.filtered[classif.filtered$transcript_id %in% c(fsm, noveltx), "FL"],
            classif.filtered[classif.filtered$transcript_id %in% novelGtx, "FL"],
            alternative = "greater")

################################################################################
# Figure S4B: Transcript Length by Category
################################################################################

len.dat <- data.frame(length = c(classif.filtered[classif.filtered$transcript_id %in% fsm, "length"],
                                 classif.filtered[classif.filtered$transcript_id %in% noveltx, "length"],
                                 classif.filtered[classif.filtered$transcript_id %in% c(fsm, noveltx), "length"],
                                 classif.filtered[classif.filtered$transcript_id %in% novelGtx, "length"]),
                      set = c(rep.int("Known", length(classif.filtered[classif.filtered$transcript_id %in% fsm, "length"])),
                              rep.int("Novel", length(classif.filtered[classif.filtered$transcript_id %in% noveltx, "length"])),
                              rep.int("Known", length(classif.filtered[classif.filtered$transcript_id %in% c(fsm, noveltx), "length"])),
                              rep.int("Novel", length(classif.filtered[classif.filtered$transcript_id %in% novelGtx, "length"]))),
                      comparison = c(rep.int("Within known genes", length(classif.filtered[classif.filtered$transcript_id %in% fsm, "length"]) + length(classif.filtered[classif.filtered$transcript_id %in% noveltx, "length"])),
                                     rep.int("Known vs novel genes", length(classif.filtered[classif.filtered$transcript_id %in% c(fsm, noveltx), "length"]) + length(classif.filtered[classif.filtered$transcript_id %in% novelGtx, "length"]))))

figs4b <- ggplot(len.dat, aes(x = comparison, y = log10(length), color = set)) +
  geom_boxplot() + 
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Transcript length (log10 bp)") +
  scale_color_manual(values = c("black", "grey"),
                     name = "Transcripts")

figs4b
ggsave("figs4b.png", figs4b, width = 4, height = 4, units = "in", dpi = 300)

# Statistical testing
wilcox.test(classif.filtered[classif.filtered$transcript_id %in% fsm, "length"],
            classif.filtered[classif.filtered$transcript_id %in% noveltx, "length"],
            alternative = "greater")

wilcox.test(classif.filtered[classif.filtered$transcript_id %in% c(fsm, noveltx), "length"],
            classif.filtered[classif.filtered$transcript_id %in% novelGtx, "length"],
            alternative = "greater")

################################################################################
# Figure S4C: Exon Count by Category
################################################################################

exons.dat <- data.frame(length = c(classif.filtered[classif.filtered$transcript_id %in% fsm, "exons"],
                                   classif.filtered[classif.filtered$transcript_id %in% noveltx, "exons"],
                                   classif.filtered[classif.filtered$transcript_id %in% c(fsm, noveltx), "exons"],
                                   classif.filtered[classif.filtered$transcript_id %in% novelGtx, "exons"]),
                        set = c(rep.int("Known", length(classif.filtered[classif.filtered$transcript_id %in% fsm, "exons"])),
                                rep.int("Novel", length(classif.filtered[classif.filtered$transcript_id %in% noveltx, "exons"])),
                                rep.int("Known", length(classif.filtered[classif.filtered$transcript_id %in% c(fsm, noveltx), "exons"])),
                                rep.int("Novel", length(classif.filtered[classif.filtered$transcript_id %in% novelGtx, "exons"]))),
                        comparison = c(rep.int("Within known genes", length(classif.filtered[classif.filtered$transcript_id %in% fsm, "exons"]) + length(classif.filtered[classif.filtered$transcript_id %in% noveltx, "exons"])),
                                       rep.int("Known vs novel genes", length(classif.filtered[classif.filtered$transcript_id %in% c(fsm, noveltx), "exons"]) + length(classif.filtered[classif.filtered$transcript_id %in% novelGtx, "exons"]))))

figs4c <- ggplot(exons.dat, aes(x = comparison, y = log10(length), color = set)) +
  geom_boxplot() + 
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Exons per transcript") +
  scale_color_manual(values = c("black", "grey"),
                     name = "Transcripts")

figs4c
ggsave("figs4c.png", figs4c, width = 4, height = 4, units = "in", dpi = 300)

# Statistical testing
wilcox.test(classif.filtered[classif.filtered$transcript_id %in% fsm, "exons"],
            classif.filtered[classif.filtered$transcript_id %in% noveltx, "exons"],
            alternative = "greater")

wilcox.test(classif.filtered[classif.filtered$transcript_id %in% c(fsm, noveltx), "exons"],
            classif.filtered[classif.filtered$transcript_id %in% novelGtx, "exons"],
            alternative = "greater")

################################################################################
# Coding Potential Analysis
################################################################################

## Import coding potential analysis results
CPC2 <- read.table("Supplementary/CPC2.txt", header = TRUE)
names(CPC2)[1] <- "transcript_id"
names(classif.filtered)[1] <- "transcript_id"
CPC2$transcript_id <- sub("\\.[^.]*$", "", CPC2$transcript_id)
CPC2 <- merge(CPC2, classif.filtered[, c("transcript_id", "structural_category")])
CPC2 <- CPC2[!is.na(CPC2$structural_category), ]

## Import ORF prediction results
TD_ORFs <- read.table("Supplementary/longest_orfs.gff3", header = FALSE)
TD_ORFs <- TD_ORFs[TD_ORFs$V3 == "CDS", c(1, 5)]
TD_ORFs <- TD_ORFs[!duplicated(TD_ORFs), ]
names(TD_ORFs) <- c("transcript_id", "ORF_length")

TD_ORFs <- TD_ORFs %>%
  group_by(transcript_id) %>%
  slice_max(order_by = ORF_length, n = 1, with_ties = FALSE) %>%
  ungroup()

TD_ORFs <- left_join(TD_ORFs, classif.filtered[, c("transcript_id", "structural_category")])
TD_ORFs <- TD_ORFs[!is.na(TD_ORFs$structural_category), ]

################################################################################
# Figure S5A: Percentage of Isoforms with ORFs by Category
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

figs5a <- ggplot(ORF_percent, aes(x = structural_category, y = percent_T,
                                  fill = structural_category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  scale_fill_manual(values = myPalette, guide = 'none') +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.title.x = element_blank()) +
  labs(y = "% of isoforms with ORFs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

figs5a
ggsave("figs5a.png", figs5a, width = 4, height = 4, units = "in", dpi = 300)

################################################################################
# Figure S5B: ORF Length by Structural Category
################################################################################

figs5b <- ggplot(TD_ORFs, aes(x = structural_category, y = log10(ORF_length + 1), fill = structural_category)) +
  geom_boxplot(outliers = FALSE) + 
  scale_fill_manual(values = myPalette) + 
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  ylab("ORF length (log10)")

figs5b
ggsave("figs5b.png", figs5b, width = 4, height = 4, units = "in", dpi = 300)

################################################################################
# Figure S5C: Coding Probability by Structural Category
################################################################################

figs5c <- ggplot(CPC2, aes(x = structural_category, y = coding_probability, fill = structural_category)) +
  geom_boxplot(outliers = FALSE) + 
  scale_fill_manual(values = myPalette) + 
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  ylab("Coding Probability")

figs5c
ggsave("figs5c.png", figs5c, width = 4, height = 4, units = "in", dpi = 300)

################################################################################
# CDS Support Analysis
################################################################################

# Get QUERY peptide lengths
query_peptides <- readAAStringSet("Assembly/ESPRESSO_corrected_SC_filtered.faa")

# Join with BLAST results
BLASTP <- read.table("Supplementary/ESPRESSO_corrected_SC_filtered_blastp.outfmt6", header = FALSE)
BLASTP <- BLASTP[, c(1, 2, 3, 4, 7, 8, 11, 12)]
names(BLASTP) <- c("transcript_id", "protein_accession", "pident", "length", "qstart", "qend", "Evalue", "score")
BLASTP <- left_join(BLASTP, query_lengths, by = "transcript_id")
BLASTP <- left_join(BLASTP, classif.filtered[, c("transcript_id", "structural_category")])
BLASTP <- BLASTP[!is.na(BLASTP$structural_category), ]

# Determine number of cores (use all but 1)
n_cores <- detectCores() - 1

# Create cluster
cl <- makeCluster(n_cores)

# Export necessary objects to cluster
query_peptides_char <- setNames(as.character(query_peptides), names(query_peptides))

clusterExport(cl, c("query_peptides_char", "BLASTP"))

# Parallel extraction
BLASTP$peptide_seq <- parSapply(cl, 1:nrow(BLASTP), function(i) {
  tid <- BLASTP$transcript_id[i]
  qstart <- BLASTP$qstart[i]
  qend <- BLASTP$qend[i]
  seq <- query_peptides_char[grepl(paste0("^", tid), names(query_peptides_char))]
  if(length(seq) == 0) return(NA)
  substring(seq[1], qstart, qend)  # Take first match only
})

# Stop cluster
stopCluster(cl)

# Add gene information
BLASTP <- left_join(BLASTP, 
                    classif.filtered[, c("transcript_id", "associated_gene")],
                    by = "transcript_id")

# For each gene, identify which peptides are unique to single transcripts
BLASTP_unique <- BLASTP %>%
  group_by(associated_gene, peptide_seq) %>%
  mutate(n_transcripts_with_peptide = n_distinct(transcript_id)) %>%
  ungroup() %>%
  filter(n_transcripts_with_peptide == 1)  # Keep only unique peptides

# Join with CPC2 coding probability scores
BLASTP_unique <- left_join(BLASTP_unique,CPC2[,c("transcript_id","coding_probability")])
BLASTP_unique <- BLASTP_unique[BLASTP_unique$coding_probability>0.5,]
BLASTP_unique <- BLASTP_unique[!is.na(BLASTP_unique$transcript_id),]

# Select best hit per transcript (highest pident, then longest match)
BLASTP_best <- BLASTP_unique %>%
  dplyr::group_by(transcript_id) %>%
  dplyr::arrange(desc(pident),desc(qend - qstart + 1)) %>%
  dplyr::slice(1) %>%
  ungroup()

# Get transcripts with unique peptide support, >50% coding probability, >95% percent identity
CDS_support <- as.character(unique(BLASTP_unique[BLASTP_unique$pident>95,"transcript_id"])$transcript_id)
BLASTP_best_pmsms <- BLASTP_best

classif.filtered$CDS_support <- FALSE
classif.filtered[classif.filtered$transcript_id %in% CDS_support, "CDS_support"] <- TRUE

classif.filtered$TD_ORF <- FALSE
classif.filtered[classif.filtered$transcript_id %in% TD_ORFs$transcript_id, "TD_ORF"] <- TRUE

CStab <- data.frame(table(classif.filtered[, c("structural_category", "TD_ORF", "CDS_support")]))
CStab <- CStab %>%
  group_by(structural_category, TD_ORF, CDS_support) %>%
  pivot_wider(names_from = CDS_support, values_from = Freq, values_fill = 0)
CStab <- data.frame(CStab)
names(CStab) <- c("structural_category", "ORF", "False", "True")
CStab$percentT <- CStab$True / (CStab$False + CStab$True)
CStab.pmsms <- CStab[10:18, ]

## BLASTP analysis with full UniProt database
BLASTP.full <- read.table("Supplementary/blastp.outfmt6", header = FALSE)
BLASTP.full <- BLASTP.full[, c(1, 2, 3, 4, 7, 8, 11, 12)]
names(BLASTP.full) <- c("transcript_id", "protein_accession", "pident", "length", "qstart", "qend", "Evalue", "score")
BLASTP.full$protein_accession <- sub("^[^|]*\\|([^|]*)\\|.*$", "\\1", BLASTP.full$protein_accession)
BLASTP.full$transcript_id <- sub("\\.[^.]*$", "", BLASTP.full$transcript_id)
BLASTP.full <- left_join(BLASTP.full, classif.filtered[, c("transcript_id", "structural_category")])
BLASTP.full <- BLASTP.full[!is.na(BLASTP.full$structural_category), ]
BLASTP.full$pident <- as.numeric(BLASTP.full$pident)

# Create cluster for parallel processing
cl <- makeCluster(n_cores)
clusterExport(cl, c("query_peptides_char", "BLASTP.full"))

# Extract peptide sequences
BLASTP.full$peptide_seq <- parSapply(cl, 1:nrow(BLASTP.full), function(i) {
  tid <- BLASTP.full$transcript_id[i]
  qstart <- BLASTP.full$qstart[i]
  qend <- BLASTP.full$qend[i]
  seq <- query_peptides_char[grepl(paste0("^", tid), names(query_peptides_char))]
  if(length(seq) == 0) return(NA)
  substring(seq[1], qstart, qend)
})
stopCluster(cl)

# Add gene information
BLASTP.full <- left_join(BLASTP.full, 
                         classif.filtered[, c("transcript_id", "associated_gene")],
                         by = "transcript_id")

# Identify unique peptides per gene
BLASTP.full_unique <- BLASTP.full %>%
  group_by(associated_gene, peptide_seq) %>%
  mutate(n_transcripts_with_peptide = n_distinct(transcript_id)) %>%
  ungroup() %>%
  filter(n_transcripts_with_peptide == 1)

# Join with CPC2 coding probability scores
BLASTP.full_unique <- left_join(BLASTP.full_unique, CPC2[, c("transcript_id", "coding_probability")])
BLASTP.full_unique <- BLASTP.full_unique[BLASTP.full_unique$coding_probability > 0.5, ]
BLASTP.full_unique <- BLASTP.full_unique[!is.na(BLASTP.full_unique$transcript_id), ]

# Select best hit per transcript (highest pident, then longest match)
BLASTP.full_unique$qstart <- as.numeric(BLASTP.full_unique$qstart)
BLASTP.full_unique$qend <- as.numeric(BLASTP.full_unique$qend)
BLASTP.full_best <- BLASTP.full_unique %>%
  dplyr::group_by(transcript_id) %>%
  dplyr::arrange(desc(pident), desc(qend - qstart + 1)) %>%
  dplyr::slice(1) %>%
  ungroup()

# Get transcripts with unique peptide support, >50% coding probability, >95% percent identity
CDS_support_full <- as.character(unique(BLASTP.full_unique[BLASTP.full_unique$pident > 95, "transcript_id"])$transcript_id)

classif.filtered$CDS_support_full <- FALSE
classif.filtered[classif.filtered$transcript_id %in% CDS_support_full, "CDS_support_full"] <- TRUE

CStab <- data.frame(table(classif.filtered[, c("structural_category", "TD_ORF", "CDS_support_full")]))
CStab <- CStab %>%
  group_by(structural_category, TD_ORF, CDS_support_full) %>%
  pivot_wider(names_from = CDS_support_full, values_from = Freq, values_fill = 0)
CStab <- data.frame(CStab)
names(CStab) <- c("structural_category", "ORF", "False", "True")
CStab$percentT <- CStab$True / (CStab$False + CStab$True)
CStab.full <- CStab[10:18, ]

################################################################################
# Figure 2C: CDS Support from Placenta MS/MS and UniProt
################################################################################
# Merge the two datasets
CStab.combined <- rbind(
  CStab.pmsms %>% mutate(source = "Placenta MS/MS"),
  CStab.full %>% mutate(source = "Full UniProt")
)

fig2c <- ggplot(CStab.combined[!CStab.combined$structural_category%in%c("Intergenic","Genic Intron"),], 
                aes(x = structural_category, y = percentT, fill = structural_category, alpha = source)) +
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_manual(values = myPalette) +
  scale_alpha_manual(values = c("Placenta MS/MS" = 1, "Full UniProt" = 0.5)) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        legend.position = "top") +
  labs(title = "Coding isoforms with unique peptides",
       subtitle = "> 50% coding probability",
       y = "% of coding isoforms") +
  scale_y_continuous(labels = scales::percent) +
  guides(fill = "none", alpha = guide_legend(override.aes = list(fill = "gray")))

fig2c
ggsave("Figures/fig2c.pdf", fig2c, width = 4.86, height = 4.18, units = "in", dpi=300)

################################################################################
# Gene Enrichment Analysis
################################################################################

## GO enrichment of top 50 genes with most isoforms
citab <- data.frame(table(classif.filtered$associated_gene))
names(citab)[1] <- "gene"
citab <- citab[rev(order(citab$Freq)), ]

tx2g.assembly <- classif.filtered[, c("transcript_id", "associated_gene")]
tx2g.assembly$transcript_id <- sub("\\.[^.]*$", "", tx2g.assembly$transcript_id)

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
# Figure S4D: Transcripts per Gene vs Exons (All Transcripts)
################################################################################

t2 <- cor.test(x = citab$isoforms, y = citab$exons)
r_val <- round(t2$estimate, 3)
p_val <- signif(t2$p.value, 3)
label_text <- paste0("r = ", r_val, "\np = ", p_val)

exoncor.df <- data.frame(isoforms = citab$isoforms,
                         exons = citab$exons)

figs4d <- ggplot(exoncor.df, aes(x = isoforms, y = exons)) +
  geom_point() +  
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal() +
  xlab("Transcripts per gene") +
  ylab("Exons per transcript") +
  annotate("text", x = 30, y = 500, label = label_text, hjust = 0) +
  ggtitle("All transcripts")

figs4d
ggsave("figs4d.png", figs4d, width = 4, height = 4, units = "in", dpi = 300)

################################################################################
# Figure S4E: Transcripts per Gene vs Exons (Highly-Expressed Transcripts)
################################################################################

t4 <- cor.test(x = citab[citab$log10TPM > 2.5, ]$isoforms, y = citab[citab$log10TPM > 2.5, ]$exons)
r_val <- round(t4$estimate, 3)
p_val <- signif(t4$p.value, 3)
label_text <- paste0("r = ", r_val, "\np = ", p_val)

exoncor.df <- data.frame(isoforms = citab[citab$log10TPM > 2.5, ]$isoforms,
                         exons = citab[citab$log10TPM > 2.5, ]$exons)

figs4e <- ggplot(exoncor.df, aes(x = isoforms, y = exons)) +
  geom_point() +  
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal() +
  xlab("Transcripts per gene") +
  ylab("Exons per transcript") +
  annotate("text", x = 30, y = 500, label = label_text, hjust = 0) +
  ggtitle("log10 TPM > 2.5")

figs4e
ggsave("figs4e.png", figs4e, width = 4, height = 4, units = "in", dpi = 300)


################################################################################
# End of Quality Control and Filtering Pipeline
################################################################################