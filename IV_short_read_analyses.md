---
output:
  html_document: default
  pdf_document: default
---

# Long-read transcriptome assembly reveals vast isoform diversity in the placenta associated with metabolic and endocrine function

## IV. Analyses with short-read datasets

**Overview**: This section presents a comprehensive analysis of short-read RNA-seq data from two independent birth cohorts (GUSTO and Gen3G) to validate and explore the biological relevance of a long-read-derived placental transcriptome assembly. Analyses include transcript quantification, global expression profiling, differential expression by gestational diabetes status (GDM), and mediation modeling of birth weight outcomes. All comparisons are made in parallel between the custom ESPRESSO-filtered long-read assembly and the GENCODE v45 reference.

**Key objectives**:

- **Quantification and validation**: Compare expression quantification accuracy and transcript detection between the ESPRESSO assembly and GENCODE v45 using dual alignment and quantification approaches.
- **Global transcriptome assessment**: Characterize transcript expression distributions, inferential uncertainty, structural features, coding potential, and inter-sample correlations across annotations and cohorts.
- **Differential expression analysis**: Identify gene-level and transcript-level associations with gestational diabetes mellitus (GDM), including cross-annotation and cross-cohort comparisons.
- **Mediation and pathway analysis**: Investigate whether specific transcripts mediate the effect of GDM on birth weight, and perform GO enrichment analysis on significant mediators to uncover functional pathways.
- **Annotation impact assessment**: Quantify the effects of transcriptome annotation choice on all analyses through parallel workflows using ESPRESSO and GENCODE references.

### A. Dataset preparation

**Purpose**: Prepare short-read RNA-seq datasets from two independent cohorts (GUSTO N=200, Gen3G N=152) for quantification against both the filtered placental assembly and GENCODE v45 reference. This dual-quantification approach enables direct comparison of transcript detection sensitivity and quantification accuracy.

### Setup

```bash
################################################################################
# Directory Structure and File Paths
################################################################################

# Primary data directories
DIR_RAW=/path/to/raw_reads                    # Raw sequencing data
DIR_TRIM=/path/to/trimmed_reads               # Quality-trimmed reads
DIR_ALIGN=/path/to/alignment_files            # BWA alignment outputs
DIR_COUNTS_ASSEMBLY=/path/to/assembly_counts  # Salmon quantification (assembly)
DIR_COUNTS_GENCODE=/path/to/gencode_counts    # Salmon quantification (GENCODE)

# Processing and temporary directories
DIR_INDEX=/path/to/preprocessing_indices_directory    # Reference genome indices
DIR_REPORTS=/path/to/reports_directory                # QC and processing reports
DIR_METRICS=/path/to/metrics_directory                # Alignment and trimming metrics

# Reference genome and annotation files
GENOME_FASTA=/path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna  # hg38 genome
REFERENCE_GTF=/path/to/gencode_v45_annotation.gtf                      # GENCODE v45

# Pre-built reference indices
GENOME_INDEX_STAR=${DIR_INDEX}/star-2.7.4a_GCA_000001405.15_GRCh38_no_alt_analysis_set
```

#### Generate Salmon mapping index for lr-assembly

**Purpose**: Create decoy-aware Salmon index for the filtered placental transcriptome assembly. Decoy sequences from the genome help distinguish between transcriptomic and genomic origins of reads, improving quantification accuracy and reducing false positive assignments.

```bash
################################################################################
# Directory Setup and Configuration
################################################################################

# Reference genome directory
DIR_GENOME=/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/GenomicReferences/genome

# Filtered assembly directory (output from SQANTI3 filtering)
DIR_TXOME=/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/STB/SQANTI3/ESPRESSO_corrected_SC_filtered

# Output directory for Salmon index
DIR_OUT=/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/STB/SQANTI3/ESPRESSO_corrected_SC_filtered/salmon_index
mkdir -p ${DIR_OUT}

# Working directory for script execution
DIR_SCRIPTS=/rsrch5/scratch/epi/stbresnahan/Placenta_Manuscript
cd ${DIR_SCRIPTS}

################################################################################
# Generate Decoy-Aware Transcriptome (Gentrome)
################################################################################

# Purpose: Create gentrome (transcriptome + genomic decoys) to improve mapping accuracy
# Decoy sequences help identify reads that might map to genomic regions not 
# represented in the transcriptome, reducing false positive transcript quantification

# Requirements:
#   - generateDecoyTranscriptome.sh (SalmonTools)
#   - bedtools
#   - salmon v1.10.2

# Activate appropriate software environment
source /rsrch5/home/epi/bhattacharya_lab/software/mashmap/bin/activate

# Generate gentrome using SalmonTools generateDecoyTranscriptome.sh
# This script extracts genomic sequences that could act as decoys for transcript mapping
sh /rsrch5/home/epi/stbresnahan/bhattacharya_lab/software/mashmap/bin/generateDecoyTranscriptome.sh \
 -j 10 \
 -b /rsrch5/home/epi/stbresnahan/bhattacharya_lab/software/bedtools \
 -a ${DIR_TXOME}/ESPRESSO_corrected_SC_filtered.gtf \
 -g /rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/STB/SQANTI3/GCA_000001405.15_GRCh38_no_alt_analysis_set_cleaned_ready.fasta \
 -t ${DIR_TXOME}/ESPRESSO_corrected_SC_filtered.fasta \
 -o ${DIR_OUT}

# Parameter explanations:
# -j 10: Use 10 threads for parallel processing
# -b: Path to bedtools binary for genomic interval operations
# -a: Assembly GTF file with transcript coordinates  
# -g: Reference genome FASTA file
# -t: Transcript sequences FASTA file
# -o: Output directory for gentrome files

# Deactivate software environment
source /rsrch5/home/epi/bhattacharya_lab/software/mashmap/bin/deactivate

################################################################################
# Build Salmon Index
################################################################################

# Purpose: Create Salmon k-mer index for fast and accurate transcript quantification
# Uses decoy-aware indexing to account for genomic sequences that could 
# interfere with transcript mapping

# Activate Salmon environment
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate salmon-1.10.2

# Build Salmon index with decoy awareness
salmon index -t ${DIR_OUT}/gentrome.fa \
  -i ${DIR_OUT}/ESPRESSO_corrected_SC_filtered \
  --decoys ${DIR_OUT}/decoys.txt \
  -k 31
  
################################################################################
# Technical Notes
################################################################################

# Parameter explanations:
# -t: Gentrome FASTA file (transcripts + decoy sequences)
# -i: Output index directory name
# --decoys: Text file listing decoy sequence names
# -k 31: k-mer size (31 is optimal for most RNA-seq applications)

# Decoy-aware indexing benefits:
# 1. Reduces false positive mappings to transcripts
# 2. Accounts for genomic regions not represented in transcriptome
# 3. Improves quantification accuracy, especially for lowly expressed transcripts
# 4. Handles reads from unannotated genomic regions

# k-mer size considerations:
# k=31 provides good balance between:
# - Specificity (longer k-mers are more specific)
# - Sensitivity (shorter k-mers capture more mappings)
# - Index size (longer k-mers create larger indices)

# Expected index characteristics:
# - Size: ~2-4 GB for human transcriptome + decoys
# - Build time: ~30-60 minutes depending on hardware
# - Memory usage: ~8-16 GB during construction
```

#### Generate Salmon mapping index for GENCODEv45

**Purpose**: Create equivalent Salmon index for GENCODE v45 reference to enable direct comparison with assembly quantification. Uses identical parameters and decoy strategy for fair comparison.

```bash
################################################################################
# Directory Setup and Configuration
################################################################################

# GENCODE reference directory
DIR_TXOME=/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/GenomicReferences/txome/gencode_v45

# Output directory for Salmon index
DIR_OUT=/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/GenomicReferences/txome/gencode.v45.salmon_index
mkdir -p ${DIR_OUT}

# Working directory for script execution
DIR_SCRIPTS=/rsrch5/scratch/epi/stbresnahan/Placenta_Manuscript
cd ${DIR_SCRIPTS}

################################################################################
# Generate Decoy-Aware Transcriptome (Gentrome) for GENCODE v45
################################################################################

# Activate appropriate software environment
source /rsrch5/home/epi/bhattacharya_lab/software/mashmap/bin/activate

# Generate gentrome using identical parameters as assembly index
sh /rsrch5/home/epi/stbresnahan/bhattacharya_lab/software/mashmap/bin/generateDecoyTranscriptome.sh \
 -j 10 \
 -b /rsrch5/home/epi/stbresnahan/bhattacharya_lab/software/bedtools \
 -a ${DIR_TXOME}/gencode.v45.annotation.gtf \
 -g /rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/STB/SQANTI3/GCA_000001405.15_GRCh38_no_alt_analysis_set_cleaned_ready.fasta \
 -t ${DIR_TXOME}/gencode.v45.transcripts.fa \
 -o ${DIR_OUT}

# Deactivate software environment
source /rsrch5/home/epi/bhattacharya_lab/software/mashmap/bin/deactivate

################################################################################
# Build Salmon Index for GENCODE v45
################################################################################

# Activate Salmon environment
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate salmon-1.10.2

# Build Salmon index with identical parameters as assembly
salmon index -t ${DIR_OUT}/gentrome.fa \
  -i ${DIR_OUT}/gencode_v45 \
  --decoys ${DIR_OUT}/decoys.txt \
  -k 31

# Consistent indexing ensures fair comparison between assembly and reference
```

#### Process GUSTO & Gen3G short reads

**Purpose**: Process short-read RNA-seq data from two independent cohorts for dual quantification against both assembly and reference transcriptomes. This provides validation across multiple datasets and population backgrounds.

**Cohort characteristics**:

- **GUSTO**: Singaporean birth cohort, term placental samples (N=200)
- **Gen3G**: French-Canadian birth cohort, term placental samples (N=152)

```bash
################################################################################
# Sample Processing Configuration
################################################################################

# Define sample identifiers for batch processing
LIBIDs=() # SampleIDs as in .fastq file names, e.g. Sample01.fastq.gz

# Processing parameters
THREADS=10  # Parallel processing threads
```

##### Adapter trimming with fastp

**Purpose**: Remove sequencing adapters and low-quality bases to improve mapping accuracy and reduce technical artifacts in quantification.

For each \$LIBID *n* in \$\{LIBIDs[1...*n*]\}
```bash
# Quality control and adapter trimming with fastp
# More stringent parameters than typical to ensure high-quality reads
fastp \
  -i ${DIR_RAW}/${LIBID}_1.fq.gz \
  -o ${DIR_TRIM}/${LIBID}_1.fastq.gz \
  -I ${DIR_RAW}/${LIBID}_2.fq.gz \
  -O ${DIR_TRIM}/${LIBID}_2.fastq.gz \
  -w ${THREADS} \
  --cut_tail \
  --cut_window_size 4 \
  --cut_mean_quality 20 \
  --length_required 36 \
  --html ${DIR_REPORTS}/${LIBID}_fastp.html \
  --json ${DIR_REPORTS}/${LIBID}_fastp.json

# Parameter explanations:
# --cut_tail: Remove low quality bases from 3' end
# --cut_window_size 4: Sliding window size for quality assessment
# --cut_mean_quality 20: Minimum average quality in sliding window
# --length_required 36: Minimum read length after trimming
# --html/--json: Generate detailed QC reports
```

##### Quantification against lr-assembly

**Purpose**: Quantify transcript expression using the filtered placental assembly. Includes inferential replicates for uncertainty assessment and bias correction for improved accuracy.

For each \$LIBID *n* in \$\{LIBIDs[1...*n*]\}
```bash
# Salmon quantification against placental assembly
TXINDEX_salmon=/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/STB/SQANTI3/ESPRESSO_corrected_SC_filtered/salmon_index/ESPRESSO_corrected_SC_filtered
DIR_COUNTS=/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/STB/infrv_GUSTO/COUNTS_assembly_SC

salmon quant \
  -i ${TXINDEX_salmon} \
  --libType A \
  --validateMappings \
  --seqBias \
  -p ${THREADS} \
  --numBootstraps 50 \
  --dumpEqWeights \
  -1 ${DIR_TRIM}/${LIBID}_1.fastq.gz \
  -2 ${DIR_TRIM}/${LIBID}_2.fastq.gz \
  -o ${DIR_COUNTS}/${LIBID}

# Parameter explanations:
# --libType A: Automatic library type detection
# --validateMappings: Use more sensitive mapping validation
# --seqBias: Correct for sequence-specific bias
# --numBootstraps 50: Generate inferential replicates for uncertainty
# --dumpEqWeights: Output equivalence class weights for analysis
```

##### Quantification against GENCODEv45

**Purpose**: Quantify the same samples against GENCODE v45 reference using identical parameters for direct comparison with assembly results.

For each \$LIBID *n* in \$\{LIBIDs[1...*n*]\}
```bash
# Salmon quantification against GENCODE v45 reference
TXINDEX_salmon=/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/GenomicReferences/txome/gencode.v45.salmon_index/gencode_v45
DIR_COUNTS=/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/STB/infrv_GUSTO/COUNTS_gencode

salmon quant \
  -i ${TXINDEX_salmon} \
  --libType A \
  --validateMappings \
  --seqBias \
  -p ${THREADS} \
  --numBootstraps 50 \
  --dumpEqWeights \
  -1 ${DIR_TRIM}/${LIBID}_1.fastq.gz \
  -2 ${DIR_TRIM}/${LIBID}_2.fastq.gz \
  -o ${DIR_COUNTS}/${LIBID}

# Identical parameters ensure fair comparison between annotations
```

#### Map GUSTO short reads for read transition analysis

**Purpose**: Perform direct read-level mapping to both transcriptomes to analyze how reads are assigned differently between assembly and reference. This provides insights into the practical impact of assembly-specific transcripts.

##### Make alignment indices for bwa-mem2

**Purpose**: Create BWA-MEM2 indices for direct read mapping analysis. BWA-MEM2 provides faster and more memory-efficient alignment than BWA-MEM while maintaining identical results.

```bash
################################################################################
# Transcriptome References for Direct Mapping
################################################################################

# Assembly transcript sequences
txome1=/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/STB/SQANTI3/ESPRESSO_corrected_SC_filtered/ESPRESSO_corrected_SC_filtered.fasta

# GENCODE v45 transcript sequences  
txome2=/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/GenomicReferences/txome/gencode_v45/gencode.v45.transcripts.fa

# Output directory for BWA indices
dir_out=/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/STB/SQANTI3/orthogonal/short-read/index

mkdir -p ${dir_out}
cd ${dir_out}

# Build BWA-MEM2 indices for both transcriptomes
bwa-mem2 index -p ESPRESSO_corrected_SC_filtered ${txome1}
bwa-mem2 index -p gencode_v45 ${txome2}
```

##### Map to transcripts with bwa-mem2

**Purpose**: Perform direct read mapping to extract read assignments for transition analysis. Filters ensure only high-quality, primary alignments are retained for analysis.

For each \$LIBID *n* in \$\{LIBIDs[1...*n*]\}
```bash
################################################################################
# Direct Read Mapping for Transition Analysis
################################################################################

# Mapping parameters optimized for transcriptome alignment
TXINDEX_BWA=/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/STB/SQANTI3/orthogonal/short-read/index
DIR_ALIGN=/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/STB/SQANTI3/orthogonal/short-read/alignments

# Map to assembly transcripts
bwa-mem2 mem -O 12,12 -E 4,4 -t ${THREADS} \
  ${TXINDEX_BWA}/ESPRESSO_corrected_SC_filtered \
  ${DIR_TRIM}/${LIBID}_1.fastq.gz ${DIR_TRIM}/${LIBID}_2.fastq.gz \
  | grep -v "^@" \
  | perl -lane 'print "$F[0]\t$F[2]" if ($F[1] & 64) and !($F[1] & 2048) and !($F[1] & 256)' \
  | gzip > ${DIR_ALIGN}/${LIBID}_assembly.tsv.gz
  
# Map to GENCODE transcripts  
bwa-mem2 mem -O 12,12 -E 4,4 -t ${THREADS} \
  ${TXINDEX_BWA}/gencode_v45 \
  ${DIR_TRIM}/${LIBID}_1.fastq.gz ${DIR_TRIM}/${LIBID}_2.fastq.gz \
  | grep -v "^@" \
  | perl -lane 'print "$F[0]\t$F[2]" if ($F[1] & 64) and !($F[1] & 2048) and !($F[1] & 256)' \
  | gzip > ${DIR_ALIGN}/${LIBID}_gencode.tsv.gz

# Parameter explanations:
# -O 12,12: Gap open penalties (higher = fewer gaps)
# -E 4,4: Gap extension penalties  
# ($F[1] & 64): First read in pair
# !($F[1] & 2048): Not supplementary alignment
# !($F[1] & 256): Not secondary alignment

# Output format: read_name \t transcript_id
# This enables downstream analysis of read assignment differences
```

### B. Global transcriptomics analyses
**Purpose**: Comprehensive comparison of quantification results between assembly and reference annotations across both cohorts. Analyzes expression patterns, statistical properties, inferential uncertainty, transcript features, coding potential, and sample correlations to validate assembly quality and biological relevance.

**Analyses performed**:

**Expression comparison**:

- **Total expression**: Sum of all transcript expression per sample
- **Mean/variance patterns**: Expression statistics by transcript categories
- **Category-specific analysis**: Shared, novel, filtered, and absent transcripts

**Statistical validation**:

- **Mixed-effects modeling**: Account for paired sample structure
- **Cross-cohort validation**: Consistent patterns across GUSTO and Gen3G
- **Uncertainty quantification**: Inferential relative variance (InfRV) analysis

**Transcript feature analysis**:

- **Full-length read support**: FL read counts by transcript categories
- **Structural characteristics**: Transcript length and exon number distributions
- **Feature correlations**: Relationship between transcript properties and expression

**Coding potential assessment**:

- **CPC2 predictions**: Coding probability scores across transcript categories
- **ORF analysis**: Open reading frame length and distribution patterns
- **Proteomics validation**: CDS support from placental mass spectrometry data
- **BLASTP homology**: Sequence similarity to known protein databases

**Quality assessment**:

- **Sample correlations**: Inter-sample consistency by transcript category
- **Technical validation**: QU correction for overdispersion
- **Biological interpretation**: Tissue-specific expression patterns

**Output figures**:

- **Main manuscript**: 1G, 1H, 1I, 2C
  * 1G: Total expression comparison across cohorts and annotations
  * 1H: Mean expression distribution by transcript categories (GUSTO)
  * 1I: Inferential uncertainty by structural category
  * 2C: Coding isoforms with CDS support from placenta MS/MS
  
- **Supplementary**: S4A-S4G, S5A-S5E, S6A-S6D
  * S4A: PCA of GUSTO samples for quality control
  * S4B: Mean/variance expression patterns (Gen3G)  
  * S4C-S4E: Inferential uncertainty analysis (both cohorts)
  * S4F-S4G: Pairwise sample correlations by transcript category
  * S5A-S5E: Transcript feature analysis (FL reads, length, exons, correlations)
  * S6A-S6D: Coding potential and ORF analysis

To generate Figures 1G-I, 2C, and S4A-G, S5A-E, S6A-D:

📜 [IVB_short-read_txomics.R](IVB_short-read_txomics.R)

### C. Differential expression of GDM
**Purpose**: Comprehensive differential expression analysis comparing gestational diabetes mellitus (GDM) cases versus controls using both gene-level (DEG) and transcript-level (DTE) approaches. Analyzes expression patterns using filtered placental transcriptome assembly versus GENCODE v45 reference across two independent cohorts.

**Analyses performed**:

**Differential expression testing**:

- **Gene-level analysis**: DESeq2 differential gene expression (DGE)
- **Transcript-level analysis**: Fishpond differential transcript expression (DTE)
- **Cross-annotation comparison**: Assembly versus GENCODE reference results
- **Cross-cohort validation**: Consistent patterns across GUSTO and Gen3G

**Statistical methodology**:

- **Batch correction**: RUVSeq with empirical control genes
- **Technical variance**: Fishpond QU (quasi-UMI) correction
- **Covariate modeling**: Mixed-effects models accounting for sex, gestational age, and batch
- **Gene-level aggregation**: Isotwas for transcript-to-gene summarization

**Case studies**:

- **CSH1 locus**: Full-splice match comparison and read mapping transitions
- **GNAS locus**: Transcript-level expression patterns by GDM status
- **Alluvial analysis**: Read mapping flow between annotations

**Cross-cohort integration**:

- **Venn diagram analysis**: Overlap of significant genes/transcripts between cohorts
- **Effect size correlation**: Log fold change comparisons across studies
- **Validation assessment**: Replication of findings between independent datasets

**Output data**:

- **GUSTO results**: GUSTO_ESPRESSO_assembly_filtered_SC.RData
- **Gen3G results**: Gen3G_ESPRESSO_assembly_filtered_SC.RData

**Output figures**:

- **Main manuscript**: 3A, 3B, 3C, 3D, 3E, 3F, 3G
  * 3A: Venn diagrams comparing DTE/DGE between cohorts and annotations
  * 3B: Log fold change comparison for CSH1 full-splice matches
  * 3C: Read count correlation for CSH1 transcripts by GDM status
  * 3D: Read mapping transitions for CSH1 between annotations
  * 3E: Log fold change comparison for GNAS full-splice matches
  * 3F: Read count correlation for GNAS transcripts by GDM status
  * 3G: Read mapping transitions for GNAS between annotations

To generate Figures 3A-G:

📜 [IVC_short-read_GDM-DE.R](IVC_short-read_GDM-DE.R)

### D. Mediation and GO enrichment analyses  
**Purpose**: Assess the potential mediating role of transcript expression in the association between gestational diabetes mellitus (GDM) and infant birth weight. Employs both single-transcript and principal component-based mediation approaches. Also performs Gene Ontology (GO) enrichment analysis of transcripts implicated in mediation using a long-read transcriptome reference.

**Analyses performed**:

**Mediation analysis**:  

- **Individual transcript-level mediation**: Using the `mediation` R package to estimate average causal mediation effects (ACME)  
- **PC-based mediation**: Dimensionality reduction on transcript expression followed by structural equation modeling (SEM) using `lavaan`  
- **Annotation-aware modeling**: Separate analyses performed using ESPRESSO long-read assembly and GENCODE v45 reference

**Gene set enrichment**:  

- **Pathway enrichment**: GO biological process enrichment via `enrichR` on transcripts with significant mediation effects  
- **Annotation comparison**: Evaluate differences in enriched pathways between ESPRESSO and GENCODE results

**Statistical methodology**:  

- **Covariate control**: Models adjust for gestational age, fetal sex, maternal ethnicity, and batch effects  
- **Significance filtering**: Multiple testing correction via Benjamini-Hochberg (FDR < 0.10 for mediation effects)  
- **Transcript filtering**: Inclusion of transcripts with sufficient expression and structural classification using SQANTI3

**Case studies**:  

- **CSH1 isoforms**: Mediation roles of individual full-splice match transcripts  
- **PC sensitivity**: Impact of PC number on mediation estimates and explained variance

**Cross-cohort integration**:  

- **GUSTO and Gen3G analyses**: Independent mediation analyses in each cohort  
- **Transcript replication**: Compare top mediating transcripts across cohorts and annotations  
- **Enrichment overlap**: Shared versus distinct GO terms across data sources

**Output data**:  

- **Individual mediation results**: `*_transcript_mediation_results.csv`  
- **PC-based SEM results**: `*_pc_mediation_results.csv`  
- **GO enrichment results**: `*_enrichR_GO_results.csv`

**Output figures**:  

- **Main manuscript**: 4B, 4C, 4D, 4E  
  * 4B: PC mediation model sensitivity analysis (GUSTO)  
  * 4C: Top transcript mediators in GUSTO (long-read assembly)  
  * 4D: Top transcript mediators in Gen3G (long-read assembly)  
  * 4E: Structural view of CSH1 transcript mediators  
- **Supplementary**: Figure S8  
  * S8: PC mediation model sensitivity analysis (Gen3G)  
- **Supplementary**: EnrichR lollipop plots for GO pathways

To generate Figures 4B–4E and S8:  

📜 [IVD_mediation_and_GOEA.R](IVD_mediation_and_GOEA.R)

---