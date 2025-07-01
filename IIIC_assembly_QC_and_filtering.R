################################################################################
# Placental Long-Read RNA-seq Assembly Quality Control, Filtering, and Characterization
# Author: Sean T. Bresnahan
# Description: Comprehensive quality control and filtering pipeline for placental
#              transcriptome assembly using SQANTI3 output. Applies multi-layered 
#              filtering based on orthogonal evidence from short-read data, CAGE-seq,
#              poly(A) sites, and artifact detection. Includes post-filtering 
#              characterization, tissue enrichment analysis, and protein validation.
# 
# Input Files:
#   - SQANTI3 classification results (ESPRESSO_corrected_classification.txt)
#   - Assembly GTF files (corrected and filtered versions)
#   - Short-read splice junction data (GUSTO_filtered_SJ.out.tab)
#   - GENCODE v45 reference annotation
#   - Assembly FASTA files (nucleotide and amino acid sequences)
#   - Junction annotations and BLASTP results
#   - Tissue-specific gene lists for enrichment analysis
#
# Output Files:
#   - Filtered classification table (PLACENTA_ASSEMBLY_classification_filtered.txt)
#   - Filtered assembly files (.gtf, .fasta, .faa, junctions)
#   - Filter reasons report (filter_reasons.txt)
#
# Output Figures:
#   - Quality Control (Pre-filtering): S2A-S2F
#     * S2A: Full-length read distribution by structural category
#     * S2B: Short-read splice junction coverage
#     * S2C: CAGE/ATAC-seq TSS validation distances
#     * S2D: Poly(A) site validation distances  
#     * S2E: Poly(A) motif distances
#     * S2F: Sequencing artifact detection summary
#   - Filtered Assembly Characterization: S3A-S3D
#     * S3A: Post-filtering transcript length distribution
#     * S3B: Post-filtering exon count distribution
#     * S3C: Known vs novel gene/isoform classification
#     * S3D: Tissue-specific enrichment of top expressed transcripts
#   - Main Manuscript Figures: 2A-2C
#     * 2A: Isoform complexity comparison (filtered assembly vs reference)
#     * 2B: CSH1 gene structure visualization
#     * 2C: Protein coding sequence support analysis
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
library(seqinr)        # Sequence analysis
library(TissueEnrich)  # Tissue-specific gene enrichment analysis
library(GSEABase)      # Gene set enrichment analysis base functions
library(data.table)    # Enhanced data manipulation
library(enrichR)       # Gene enrichment analysis interface
library(org.Hs.eg.db)  # Human genome annotation database
library(ggtranscript)  # Transcript structure visualization

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
classif <- read.table("ESPRESSO_corrected_classification.txt", header = TRUE)

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
GTF <- import("ESPRESSO_corrected_corrected.gtf")

# Import short-read splice junction data from STAR aligner
# This provides orthogonal validation of splice sites
SJ.tab <- read.table("GUSTO_filtered_SJ.out.tab", header = FALSE, sep = "\t")

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
  summarise(
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
  group_by(transcript_id) %>%
  arrange(seqnames, start, end, .by_group = TRUE) %>%
  mutate(exon_number = row_number()) %>%
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
  group_by(transcript_id, exon_number) %>%
  summarise(total_reads_sum = sum(total_reads), .groups = "drop")

# Calculate TSS ratio: first exon reads / last exon reads
# Lower ratios may indicate degradation or incomplete transcripts
tss_ratio <- summarized_exons %>%
  group_by(transcript_id) %>%
  summarise(
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
reasons <- rbind(
  # filt.cov,  # Missing in provided code
  filt.5prime,
  filt.3prime,
  filt.perc_A_downstream_TTS,
  filt.RTS_stage,
  filt.predicted_NMD
)

# Write filtering report
write.table(reasons, "/filter_reasons.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Create filtered classification table
classif.filtered <- classif[!classif$isoform %in% reasons$isoform, ]
write.table(classif.filtered, "PLACENTA_ASSEMBLY_classification_filtered.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

################################################################################
# Filter Supporting Files
################################################################################

## Filter protein sequences (amino acid FASTA)
faa <- read.fasta(file = "PLACENTA_ASSEMBLY_corrected_filtered.faa", 
                  seqtype = "AA", as.string = TRUE, set.attributes = FALSE)

# Extract transcript IDs from FASTA headers and filter
faa <- faa[sapply(strsplit(names(faa), "\t"), `[[`, 1) %in% classif.filtered$isoform]
write.fasta(sequences = faa, names = names(faa), 
            file.out = paste0(DIR_OUT, "/", OUT_PREFIX, ".faa"))
rm(faa)

## Filter nucleotide sequences (DNA FASTA)
fasta <- read.fasta(file = "PLACENTA_ASSEMBLY_corrected_filtered.fasta", 
                    seqtype = "DNA", as.string = TRUE, set.attributes = FALSE)
fasta <- fasta[names(fasta) %in% classif.filtered$isoform]
write.fasta(sequences = fasta, names = names(fasta), 
            file.out = paste0(DIR_OUT, "/", OUT_PREFIX, ".fasta"))
rm(fasta)

## Filter splice junction annotations
junc <- read.table("PLACENTA_ASSEMBLY_junctions.txt", header = TRUE)
junc <- junc[junc$isoform %in% classif.filtered$isoform, ]
write.table(junc, "PLACENTA_ASSEMBLY_junctions_filtered.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

## Filter GTF annotation
gtf <- readGFF("PLACENTA_ASSEMBLY_corrected.gtf")
gtf <- gtf[gtf$transcript_id %in% classif.filtered$isoform, ]

# Format GTF attributes for standard compliance
gtf$score <- "."
gtf$phase <- "."
gtf$info <- NA

# Create properly formatted attribute strings
gtf[gtf$type == "transcript", "info"] <- paste0("transcript_id \"",
                                                gtf[gtf$type == "transcript", "transcript_id"],
                                                "\"; gene_id \"",
                                                gtf[gtf$type == "transcript", "gene_id"], "\"")

gtf[gtf$type == "exon", "info"] <- paste0("transcript_id \"",
                                          gtf[gtf$type == "exon", "transcript_id"],
                                          "\"; gene_id \"",
                                          gtf[gtf$type == "exon", "gene_id"], "\";")

# Remove redundant columns
gtf$transcript_id <- NULL
gtf$gene_id <- NULL

# Write filtered GTF
write.table(gtf, "PLACENTA_ASSEMBLY_corrected_filtered.gtf", 
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
ggsave("figs2a.png", figs2a, width = 4, height = 4, units = "in", dpi = 300)

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
ggsave("figs2b.png", figs2b, width = 4, height = 4, units = "in", dpi = 300)

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
ggsave("figs2c.png", figs2c, width = 3, height = 3, units = "in", dpi = 300)

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
ggsave("figs2d.png", figs2d, width = 3, height = 3, units = "in", dpi = 300)

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
ggsave("figs2e.png", figs2e, width = 3, height = 3, units = "in", dpi = 300)

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
ggsave("figs2f.png", figs2f, width = 6, height = 3, units = "in", dpi = 300)

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
ggsave("figs3a.png", figs3a, width = 3, height = 3, units = "in", dpi = 300)

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
ggsave("figs3b.png", figs3b, width = 4, height = 4, units = "in", dpi = 300)

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
ggsave("figs3c.png", figs3c, width = 3, height = 3, units = "in", dpi = 300)

################################################################################
# Figure S3D: Tissue-Specific Enrichment Analysis
################################################################################

# Regenerate transcript-to-gene mapping for filtered assembly
tx2g.assembly <- data.frame(import("PLACENTA_ASSEMBLY_corrected_filtered.gtf"))
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

# Import GENCODE v45 reference for comparison
tx2g.gencode <- data.frame(import("gencode.v45.annotation.gtf"))
tx2g.gencode <- tx2g.gencode[tx2g.gencode$seqnames %in% 
                               c(paste0("chr", seq.int(1, 22, 1)), "chrX", "chrY", "chrM"), ]
tx2g.gencode <- tx2g.gencode[, c("transcript_id", "gene_id")]
tx2g.gencode <- tx2g.gencode[!is.na(tx2g.gencode$transcript_id), ]
tx2g.gencode <- tx2g.gencode[!duplicated(tx2g.gencode), ]

# Save tx2g objects for later use
save(list=c("tx2g.assembly","tx2g.gencode"),file="tx2g.RData")

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
ggsave("figs3d.png", figs3d, width = 10.55, height = 5.92, units = "in", dpi = 300)

################################################################################
# Figure 2A: Isoform Complexity Comparison (Filtered Assembly)
################################################################################

# Prepare filtered assembly transcript-to-gene mapping
tx2g.assembly.filtered <- classif.filtered[, c("transcript_id", "associated_gene")]
names(tx2g.assembly.filtered)[2] <- "gene_id"

# Calculate isoforms per gene for filtered assembly
iso.tab.assembly.filtered <- data.frame(table(tx2g.assembly.filtered$gene_id))
names(iso.tab.assembly.filtered) <- c("gene_id", "Freq.assembly")
iso.tab.assembly.filtered$gene_id <- as.character(iso.tab.assembly.filtered$gene_id)

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
  xlab("GENCODEv45 (x10³)") +
  ylab("lr-assembly (x10³)") +
  ggtitle("Isoforms per gene")

# Add marginal histograms
ggMarginal(fig2a, type = "histogram", fill = "grey40", color = "white")

################################################################################
# Figure 2B: CSH1 Gene Structure Visualization
################################################################################

# Prepare combined GTF data for transcript structure plotting
gtf.assembly <- data.frame(import("PLACENTA_ASSEMBLY_corrected_filtered.gtf"))
gtf.assembly$gene_id <- NULL
gtf.assembly <- left_join(gtf.assembly, tx2g.assembly)
gtf.assembly <- left_join(gtf.assembly, classif[, c("transcript_id", "structural_category")])
names(gtf.assembly)[12] <- "class"

# Import reference transcripts for comparison
gtf.gencode <- data.frame(import("gencode.v45.annotation.gtf"))
gtf.gencode <- gtf.gencode[, names(gtf.gencode) %in% names(gtf.assembly)]
gtf.gencode$class <- "Reference"
gtf.gencode <- gtf.gencode[!gtf.gencode$transcript_id %in% gtf.assembly$transcript_id, ]

# Combine assembly and reference annotations
gtf.assembly <- rbind(gtf.assembly, gtf.gencode)
gtf.assembly$gene_id <- sapply(strsplit(gtf.assembly$gene_id, '[.]'), function(x) x[1])

# Function to create transcript structure plots
tx.plot <- function(gtf.df, gene.name, title) {
  # Filter data for specific gene and exons only
  dat <- gtf.df[gtf.df$gene_id == gene.name & gtf.df$type == "exon", ]
  dat$transcript_id <- factor(dat$transcript_id, levels = unique(dat$transcript_id))
  
  # Define structural category hierarchy for ordering
  dat$class <- factor(as.character(dat$class),
                      levels = c("Reference", "FSM", "ISM", "NIC", "NNC",
                                 "Genic Genomic", "Antisense", "Fusion",
                                 "Intergenic", "Genic Intron"))
  
  # Order transcripts by structural category
  dat <- dat %>%
    mutate(transcript_id = factor(transcript_id,
                                  levels = dat %>%
                                    distinct(transcript_id, class) %>%
                                    arrange(class) %>%
                                    pull(transcript_id)))
  
  # Create transcript structure plot
  p <- dat %>% 
    ggplot(aes(xstart = start, xend = end, y = transcript_id)) +
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

# Create CSH1 gene structure plot
fig2b <- tx.plot(gtf.assembly, "ENSG00000136488", "CSH1 (ENSG00000136488)") +
  scale_fill_manual(values = c("#66c2a5", "#8da0cb", "#e78ac3"))

fig2b
ggsave("fig2b.png", fig2b, width = 6, height = 8, units = "in", dpi = 300)

################################################################################
# Figure 2C: Protein Coding Sequence Support Analysis
################################################################################

# Import BLASTP results for protein homology validation
blastp <- read.table("ESPRESSO_corrected_SC_filtered_blastp.outfmt6", header = FALSE, sep = "\t")

# Filter for high-confidence protein matches (≥70% identity)
blastp <- unique(blastp[blastp$V3 >= 70, 1])

# Prepare data for protein support analysis
fig2c.dat0 <- classif.filtered[, c("transcript_id", "structural_category", "coding")]

# Classify isoforms as known (FSM) or novel
fig2c.dat0$isoform_class <- "Novel"
fig2c.dat0[fig2c.dat0$transcript_id %in% classif.filtered[classif.filtered$structural_category == "FSM", "transcript_id"], "isoform_class"] <- "Known"

# Classify genes as known (ENSEMBL) or novel
fig2c.dat0$gene_class <- "Known"
fig2c.dat0[fig2c.dat0$transcript_id %in% classif.filtered[grep("novel", classif.filtered$associated_gene), "transcript_id"], "gene_class"] <- "Novel"

# Mark transcripts with human protein database support
fig2c.dat0$human_protein_support <- FALSE
fig2c.dat0[fig2c.dat0$transcript_id %in% blastp, "human_protein_support"] <- TRUE

# Filter out inconsistent combinations for cleaner analysis
fig2c.dat0$keep <- TRUE

# Remove known coding isoforms without protein support (likely false positives)
fig2c.dat0[fig2c.dat0$coding == "coding" & fig2c.dat0$isoform_class == "Known" & fig2c.dat0$human_protein_support == FALSE, "keep"] <- FALSE

# Remove non-coding transcripts with protein support (likely false positives)  
fig2c.dat0[fig2c.dat0$coding == FALSE & fig2c.dat0$human_protein_support == TRUE, "keep"] <- FALSE

# Focus analysis on coding transcripts only
fig2c.dat0 <- fig2c.dat0[fig2c.dat0$keep == TRUE, ]
fig2c.dat0 <- fig2c.dat0[fig2c.dat0$coding == "coding", ]

# Count combinations and calculate percentages
fig2c.dat <- data.frame(table(fig2c.dat0[, c("isoform_class", "gene_class", "human_protein_support")]))
fig2c.dat$percent <- 0

# Calculate observed percentages from the data
# These values appear to be hardcoded based on actual counts
fig2c.dat[2, 5] <- 932 / (932 + 5912)     # Known isoform, known gene, with protein support
fig2c.dat[4, 5] <- 83 / (83 + 24)         # Novel isoform, known gene, with protein support  
fig2c.dat[5, 5] <- 1                      # Novel isoform, novel gene, with protein support
fig2c.dat[6, 5] <- 5912 / (932 + 5912)    # Known isoform, known gene, without protein support
fig2c.dat[8, 5] <- 24 / (83 + 24)         # Novel isoform, known gene, without protein support

# Remove empty combinations
fig2c.dat <- fig2c.dat[-c(3, 7), ]

# Create protein support visualization
fig2c <- ggplot(fig2c.dat, aes(y = percent, x = gene_class, fill = isoform_class)) +
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_manual(values = c("black", "grey")) +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "top") +
  labs(title = "Known vs novel with\nputative proteins",
       y = "% of coding isoforms", 
       x = "Gene",
       fill = "Isoform")

fig2c
ggsave("fig2c.png", fig2c, width = 3, height = 3, units = "in", dpi = 300)

################################################################################
# End of Quality Control and Filtering Pipeline
################################################################################