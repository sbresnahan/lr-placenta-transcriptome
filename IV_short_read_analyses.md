---
output:
  html_document: default
  pdf_document: default
---

# Long-read assembly reveals vast transcriptional complexity in the placenta associated with metabolic and endocrine function

## IV. Analyses with short-read datasets

**Overview**: This section presents a comprehensive analysis of short-read RNA-seq data from two independent birth cohorts (GUSTO and Gen3G) to validate and explore the biological relevance of a long-read-derived placental transcriptome assembly. Analyses include transcript quantification, global expression profiling, differential expression by gestational diabetes status (GDM), and mediation modeling of birth weight outcomes. All comparisons are made in parallel between the custom ESPRESSO-filtered long-read assembly and the GENCODE v45 reference.

*Note: several analyses performed here were repeated for the GTEx v9 long-read and short-read data,
which are available under controlled access via dbGaP [accession phs000424.v9](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000424.v9) 
and on [AnVIL](https://anvil.terra.bio/#workspaces/anvil-datastorage/AnVIL_GTEx_V9_hg38).
See the manuscript and [orthogonal_for_GTEx.txt](Datasets/orthogonal_for_GTEx.txt) for details.*

**Objectives**:

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

#### Generate combined GTF (lr-assembly + GENCODE v45)

**Purpose**: Create a combined GTF annotation that merges novel transcripts and genes from the filtered long-read assembly with the full GENCODE v45 reference. This combined annotation (GENCODE+) enables quantification that captures both known and novel placental isoforms in a single reference, allowing direct comparison with assembly-only and GENCODE-only results.

**Process**:

1. Read the filtered assembly GTF and classification table
2. Resolve gene IDs for novel transcripts using SQANTI3 associated gene assignments
3. Create gene IDs for truly novel genes (no GENCODE association)
4. Read GENCODE v45 GTF and format to match assembly structure
5. Identify novel genes and novel transcripts not present in GENCODE
6. Merge GENCODE v45 with assembly-specific entries (novel genes and novel isoforms of known genes)
7. Sort combined GTF by genomic coordinates with transcripts preceding their exons

📜 [combine_GTF.R](combine_GTF.R)

#### Generate Salmon mapping index for combined GTF (GENCODE+)

**Purpose**: Create a decoy-aware Salmon index for the combined lr-assembly + GENCODE v45 annotation. Uses identical indexing parameters as the assembly-only and GENCODE-only indices for fair comparison across all three annotation strategies.

📜 [make_salmon_combined.lsf](make_salmon_combined.lsf)

```bash
################################################################################
# Directory Setup and Configuration
################################################################################

DIR_GENOME=/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/genome
DIR_TXOME=/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/combined_GTF
DIR_OUT=/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/combined_GTF/salmon
mkdir -p ${DIR_OUT}

################################################################################
# Generate Decoy-Aware Transcriptome and Build Salmon Index
################################################################################

# Steps (commented sections in script show full gentrome generation):
# 1. Extract transcript sequences from combined GTF using gffread
# 2. Generate decoy transcriptome using generateDecoyTranscriptome.sh
# 3. Build Salmon index with decoy awareness

eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate salmon-1.10.2

cd ${DIR_OUT}

salmon index -t gentrome.fa \
  -i combined \
  --decoys decoys.txt \
  -k 31
```

#### Quantification against combined annotation (GENCODE+)

**Purpose**: Quantify short-read RNA-seq from both cohorts against the combined lr-assembly + GENCODE v45 annotation. This provides the GENCODE+ quantification used in downstream comparisons alongside assembly-only and GENCODE-only results.

##### GUSTO quantification against combined annotation

📜 [run_salmon_combined_GUSTO.lsf](run_salmon_combined_GUSTO.lsf) | 📜 [loop_salmon_combined_GUSTO.sh](loop_salmon_combined_GUSTO.sh)

For each \$LIBID *n* in \$\{LIBIDs[1...*n*]\}
```bash
TXINDEX_salmon=/rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/combined_GTF/salmon/combined
DIR_TRIM=/rsrch5/home/epi/bhattacharya_lab/data/mapqtl/GUSTO/rnaseq/fastq
DIR_COUNTS=/rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/infrv_GUSTO/COUNTS_combined
mkdir -p ${DIR_COUNTS}

salmon quant \
  -i ${TXINDEX_salmon} --libType A \
  --validateMappings --seqBias \
  -p 12 --numBootstraps 100 \
  -1 ${DIR_TRIM}/${LIBID}_1.fq.gz \
  -2 ${DIR_TRIM}/${LIBID}_2.fq.gz \
  -o ${DIR_COUNTS}/${LIBID}
```

##### Gen3G quantification against combined annotation

📜 [run_salmon_combined_Gen3G.lsf](run_salmon_combined_Gen3G.lsf) | 📜 [loop_salmon_combined_Gen3G.sh](loop_salmon_combined_Gen3G.sh)

For each \$LIBID *n* in \$\{LIBIDs[1...*n*]\}
```bash
TXINDEX_salmon=/rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/combined_GTF/salmon/combined
DIR_TRIM=/rsrch5/home/epi/bhattacharya_lab/data/mapqtl/Gen3G
DIR_COUNTS=/rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/infrv_Gen3G/COUNTS_combined
mkdir -p ${DIR_COUNTS}

salmon quant \
  -i ${TXINDEX_salmon} --libType A \
  --validateMappings --seqBias \
  -p 12 --numBootstraps 100 \
  -1 ${DIR_TRIM}/${LIBID}_1.fastq \
  -2 ${DIR_TRIM}/${LIBID}_2.fastq \
  -o ${DIR_COUNTS}/${LIBID}
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
  --numBootstraps 100 \
  --dumpEqWeights \
  -1 ${DIR_TRIM}/${LIBID}_1.fastq.gz \
  -2 ${DIR_TRIM}/${LIBID}_2.fastq.gz \
  -o ${DIR_COUNTS}/${LIBID}

# Parameter explanations:
# --libType A: Automatic library type detection
# --validateMappings: Use more sensitive mapping validation
# --seqBias: Correct for sequence-specific bias
# --numBootstraps 100: Generate inferential replicates for uncertainty
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
  --numBootstraps 100 \
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

### B. Short-Read Transcriptomics Validation

**Purpose**: Comprehensive validation of the filtered placental transcriptome assembly versus GENCODE v45 reference using short-read RNA-seq quantification from two independent cohorts (GUSTO n=200, Gen3G n=152). Analyzes expression patterns, quantification uncertainty, transcript features, and coding potential to demonstrate assembly quality and biological utility.

**Cohorts**:
- **GUSTO (Singapore)**: N=200 placental samples, fetal-facing biopsies
- **Gen3G (Canada)**: N=152 placental samples, fetal-facing biopsies
- **Quality control**: PCA-based outlier removal, batch effect filtering, RIN score assessment

**Annotations compared**:
- **GENCODEv45**: Reference annotation only
- **GENCODE+**: Reference + filtered assembly combined
- **lr-assembly**: Filtered placental assembly only

**Quantification methods**:
- **Salmon**: Transcript-level pseudo-alignment and quantification
- **QU correction**: Quasi-UMI correction for overdispersion using fishpond's catchSalmon
- **Inferential replicates**: Bootstrap-based uncertainty estimation via tximeta
- **Normalization**: DESeq2 size factor normalization, TPM calculations

**Transcript classification**:
- **Shared**: Transcripts present in both assembly and GENCODE
- **Novel**: Assembly-specific transcripts from known genes
- **Filtered**: Transcripts removed during quality control
- **Absent**: GENCODE transcripts not assembled

**Analyses performed**:

**1. Expression comparison**:
- Total isoform expression per sample (sum of log TPM)
- Mean expression distributions by transcript category
- Variance patterns across transcript categories
- Cross-cohort consistency validation

**2. Inferential uncertainty analysis (InfRV)**:
- Quantification uncertainty by annotation (GENCODEv45 vs GENCODE+ vs lr-assembly)
- Uncertainty by structural category (FSM, ISM, NIC, NNC, Other)
- Relationship between uncertainty and exon sharing (overlap analysis)
- Cross-tissue comparison (Placenta vs GTEx Adipose vs GTEx Fibroblasts)

**3. Statistical validation**:
- Mixed-effects models accounting for paired sample structure
- Pairwise contrasts with Tukey adjustment
- Estimated marginal means comparisons
- Cross-cohort replication

**4. Sample quality and correlation**:
- Principal component analysis (PCA) for batch effects
- Pairwise sample correlations by structural category
- Transcriptional complexity by fetal sex
- Technical reproducibility assessment

**5. Exon sharing analysis**:
- Overlap counting with ≥80% reciprocal overlap threshold
- Parallel computation for efficiency
- Correlation with inferential uncertainty
- Regression statistics (R, R², p-values)

**6. Cross-tissue validation**:
- GUSTO samples quantified against GTEx tissue assemblies
- Adipose - Subcutaneous assembly
- Cells - Cultured Fibroblasts assembly
- Inferential uncertainty comparison across tissues

**Output files**:

**Data files**:
- **tx2g_ESPRESSO_assembly_SC_filtered.RData**: Transcript-to-gene mappings (assembly, gencode, combined)
- **metadata_GUSTO.RData**: GUSTO cohort metadata and covariates
- **metadata_Gen3G.RData**: Gen3G cohort metadata and covariates
- **infrv_GUSTO.RDS**: Inferential relative variance data for GUSTO

**Quality control figures**:
- **Figure S6A**: PCA plots for GUSTO (25.72% PC1, 18.26% PC2) and Gen3G (19.54% PC1, 6.36% PC2) cohorts
- **tx_complex_by_sex.png**: Transcriptional complexity (isoforms per gene) by fetal sex

**Expression comparison figures**:
- **Figure 3B**: Total isoform expression across samples (violin + boxplot, GUSTO and Gen3G, 3 annotations)
- **Figure 3C**: Mean and variance by transcript categories (GUSTO, grouped by Discovery vs Filtration)
- **Figure S6B**: Mean and variance by transcript categories (Gen3G replication)

**Inferential uncertainty figures**:
- **Figure 3D**: InfRV distribution by annotation (GUSTO: GENCODEv45, GENCODE+, lr-assembly) with mode indicators
- **Figure S6C**: InfRV distribution by annotation (Gen3G replication)
- **Figure S6D**: InfRV by structural category (GUSTO, lr-assembly only: FSM, ISM, NIC, NNC, Other)
- **Figure S6E**: InfRV by structural category (Gen3G replication)
- **Figure S6H**: InfRV vs exon sharing scatter plot with linear regression (R, R², p-values for each annotation)
- **Figure S6I**: InfRV across placenta and GTEx tissues (Placenta, Fibroblasts, Adipose vs GENCODEv45, GENCODE+)

**Sample correlation figures**:
- **Figure S6F**: Pairwise sample correlations by structural category (GUSTO, lr-assembly)
- **Figure S6G**: Pairwise sample correlations by structural category (Gen3G, lr-assembly)

To generate all short-read transcriptomics validation figures and analysis:
📜 [IVB_short-read_txomics.R](IVB_short-read_txomics.R)

### C. Differential Expression Analysis: GDM Association Study

**Purpose**: Comprehensive differential expression analysis comparing gestational diabetes mellitus (GDM) cases versus controls using both gene-level (DEG) and transcript-level (DTE) approaches across two independent cohorts (GUSTO n=200, Gen3G n=152). Compares results from filtered placental transcriptome assembly, GENCODE v45 reference, and GENCODE+ combined annotation to demonstrate the utility of tissue-specific isoform discovery.

**Cohorts and sample sizes**:
- **GUSTO (Singapore)**: N=200 placental samples (discovery cohort)
- **Gen3G (Canada)**: N=152 placental samples (validation cohort)
- **Cross-tissue analysis**: GUSTO samples quantified against GTEx Adipose and Fibroblasts assemblies

**Annotations compared**:
- **GENCODEv45**: Reference annotation only
- **GENCODE+**: Reference + filtered assembly combined
- **lr-assembly**: Filtered placental assembly only
- **GTEx tissues**: Adipose - Subcutaneous and Cells - Cultured Fibroblasts assemblies

**Statistical methodology**:

**1. Batch correction and normalization**:
- **RUVSeq (Remove Unwanted Variation)**: RUVr method with k=2 factors
- Residuals from GLM fit used to identify unwanted variation
- PCA visualization before and after batch correction
- Integration of W1 and W2 factors into differential expression model

**2. Quantification preprocessing**:
- **Salmon pseudo-alignment**: Transcript-level quantification
- **QU correction**: fishpond's catchSalmon for overdispersion adjustment
- **DESeq2 normalization**: Size factor estimation and variance stabilization
- **Gene-level aggregation**: tximport summarization from transcripts

**3. Differential expression models**:
- **Full model**: `~ sex + GA + W1 + W2 + trait` (with RUV correction)
- **Reduced model**: `~ sex + GA + trait` (without RUV)
- **Trait**: GDM status (categorical: 0 vs 1)
- **Covariates**: Fetal sex, gestational age (GA), batch factors

**4. Multiple testing correction**:
- **FDR control**: Benjamini-Hochberg adjustment
- **Significance threshold**: padj < 0.1 for DTE/DEG identification
- **Gene-level screening**: isotwas p_screen for transcript-to-gene aggregation
- **Confirmation**: isotwas p_confirm with adaptive alpha threshold

**Cross-annotation comparison strategies**:

**1. Transcript categorization**:
- **Shared**: Present in both assembly and GENCODE
- **Novel**: Assembly-specific transcripts from known genes
- **Filtered**: Transcripts removed during quality control (absent from lr-assembly, present in unfiltered)
- **Absent**: GENCODE transcripts not detected in assembly

**2. Gene set overlaps**:
- DTE/DEG overlap across cohorts (GUSTO, Gen3G)
- DTE/DEG overlap across annotations (GENCODEv45, GENCODE+, lr-assembly)
- Cross-tissue specificity (Placenta vs Adipose vs Fibroblasts)

**Case study analyses**:

**1. CSH1 (Chorionic Somatomammotropin Hormone 1)**:
- 6 full-splice match (FSM) transcripts shared between annotations
- Log fold change comparison between lr-assembly and GENCODEv45
- Read count correlation highlighting ENST00000558284.1
- Alluvial diagram tracking read mapping transitions from GENCODEv45 → lr-assembly
- Separate analysis by GDM status (cases vs controls)

**2. GNAS (GNAS Complex Locus)**:
- 77 FSM transcripts in GENCODE, subset detected in assembly
- Log fold change comparison focusing on ENST00000476196.5
- Read count correlation by GDM status
- Alluvial flow diagram showing redistribution of reads across structural categories
- Quantification of read reassignment from non-GNAS to GNAS isoforms

**Read mapping transition analysis**:
- **Parallel processing**: data.table and doParallel for efficient computation
- **Sample selection**: 5 GDM=0 and 5 GDM=1 samples (random sampling, seed=1234)
- **Read-to-transcript mapping**: Separate files for GENCODEv45 and lr-assembly
- **Transition categories**: Unmapped, target gene, non-target, structural categories
- **Frequency scaling**: Square root transformation for alluvial visualization
- **Statistical summaries**: Percentage of reads redistributed across annotations

**Visualization methods**:
- **Volcano plots**: Log fold change vs -log10(adjusted p-value) with point thinning (200x200 bins, max 5 points/bin)
- **UpSet plots**: Intersection size plots with colored combination matrix
- **Error bar plots**: Log fold changes with 95% confidence intervals (±1.96 × SE)
- **Scatter plots**: Read count correlations with GDM status highlighting
- **Alluvial diagrams**: Read mapping flow between annotations using ggalluvial

**Output files**:

**Data files**:
- **GUSTO_ESPRESSO_assembly_filtered_SC.RData**: All GUSTO differential expression results
  - Assembly, GENCODE, GENCODE+, Adipose, Fibroblasts
  - Gene-level and transcript-level results
  - VST-transformed counts for visualization
- **Gen3G_ESPRESSO_assembly_filtered_SC.RData**: All Gen3G differential expression results
  - Assembly, GENCODE, GENCODE+ annotations
  - Gene-level and transcript-level results
- **TableS10.csv**: Differential transcript expression results across all cohorts and annotations
- **TableS11.csv**: Differential gene expression results across all cohorts and annotations

**Differential expression overview figures**:
- **Figure S7A**: DTE volcano plots (6 panels: GUSTO/Gen3G × GENCODEv45/GENCODE+/lr-assembly)
- **Figure S7B**: DEG volcano plots (6 panels: GUSTO/Gen3G × GENCODEv45/GENCODE+/lr-assembly)
- **Figure 3E**: UpSet plot comparing DTE and DEG across cohorts and annotations (12 sets total, filtered to n≥10)
- **Figure S8C**: UpSet plot for GUSTO tissues (Placenta, Adipose, Fibroblasts; 10 sets total, filtered to n≥10)

**CSH1 case study figures**:
- **Figure 5B**: Log fold change comparison for CSH1 FSM transcripts (lr-assembly vs GENCODEv45)
  - Error bars showing 95% CI for 6 FSM isoforms
  - ENST00000558284.1 highlighted in black with label
  - Reference lines at logFC=0 for both axes
- **Figure S5C**: Read count correlation for CSH1 (log2 scale, 0-25 range)
  - ENST00000558284.1 colored by GDM status (red=GDM, blue=control)
  - All other CSH1 FSM in grey
  - Diagonal reference line (slope=1, intercept=0)
- **Figure 5D**: Alluvial diagram for CSH1 read mapping transitions
  - Left axis: GENCODEv45 transcripts (6 CSH1 FSM + non-CSH1)
  - Right axis: lr-assembly categories (ENST00000558284.1, FSM, ISM, NIC, NNC, non-CSH1)
  - Y-axis: Square root of read count
  - Combined data from 5 GDM=0 and 5 GDM=1 samples

**GNAS case study figures**:
- **Figure S5E**: Log fold change comparison for GNAS FSM transcripts (lr-assembly vs GENCODEv45)
  - Error bars showing 95% CI for detected FSM isoforms
  - ENST00000476196.5 highlighted in black with label
  - Reference lines at logFC=0 for both axes
- **Figure S5F**: Read count correlation for GNAS (log2 scale, 0-12 range)
  - ENST00000476196.5 colored by GDM status (red=GDM, blue=control)
  - All other GNAS FSM in grey
  - Diagonal reference line
- **Figure 5G**: Alluvial diagram for GNAS read mapping transitions
  - Left axis: GENCODEv45 transcripts (ENST00000476196.5, Other GNAS n=76, non-GNAS)
  - Right axis: lr-assembly categories (ENST00000476196.5, FSM, non-GNAS)
  - Frequency adjustment: 2× for ENST00000476196.5, 8× for non-GNAS
  - Statistical summary of read redistribution percentage

To generate all differential expression results, cross-cohort comparisons, and case study visualizations:
📜 [IVC_DGE_DTE.R](IVC_DGE_DTE.R)

### D. Mediation Analysis and Gene Ontology Enrichment

**Purpose**: Investigate how differentially expressed placental transcripts mediate the relationship between gestational diabetes mellitus (GDM) and infant birth weight. Employs both individual transcript-level and principal component (PC)-based mediation approaches to assess indirect effects. Performs Gene Ontology (GO) enrichment analysis to identify biological processes enriched among significant mediators and differentially expressed genes/transcripts.

**Cohorts and sample sizes**:
- **GUSTO (Singapore)**: N=200 placental samples
- **Gen3G (Canada)**: N=152 placental samples
- **Cross-tissue analysis**: GUSTO samples with GTEx Adipose and Fibroblasts assemblies

**Annotations compared**:
- **GENCODEv45**: Reference annotation only
- **GENCODE+**: Reference + filtered assembly combined
- **lr-assembly**: Filtered placental assembly only
- **GTEx tissues**: Adipose - Subcutaneous and Cells - Cultured Fibroblasts

**Individual transcript mediation analysis**:

**1. Statistical framework**:
- **Treatment (X)**: GDM status (0 vs 1)
- **Mediator (M)**: Variance-stabilized transcript expression (VST)
- **Outcome (Y)**: Birth weight (scaled)
- **Path a**: Effect of GDM on transcript expression (mediator model: `TXE ~ GDM + sex + GA`)
- **Path b**: Effect of transcript expression on birth weight, adjusting for GDM (outcome model: `birth_weight ~ GDM + TXE + sex + GA`)
- **Indirect effect**: a × b (average causal mediation effect, ACME)

**2. Bootstrap inference**:
- Bias-corrected confidence intervals via percentile method
- P-value threshold: p < 0.05 for significance

**3. Parallel processing**:
- Foreach + doParallel for computational efficiency
- Parallel analysis of significant differentially expressed transcripts (padj < 0.1)
- Error handling for failed mediation models

**4. Union transcript analysis**:
- Top 10 mediators from GUSTO + top 10 from Gen3G
- Re-analysis of union set (N=20 unique transcripts) in both cohorts
- Cross-cohort validation of mediation effects

**Principal component mediation analysis**:

**1. Dimensionality reduction**:
- PCA on highly significant transcripts (FDR < 0.1 from differential expression)
- Separate PC analysis for each annotation (GENCODEv45, GENCODE+, lr-assembly)
- Variance explained by each PC and cumulative variance
- PC-ethnicity correlation analysis (ANOVA with eta-squared effect size)

**2. Structural equation modeling (SEM)**:
- lavaan package for path analysis
- **Mediator equations**: `PC_i ~ a_i*GDM + sex_i*sex + GA_i*GA + ethnicity_i*ethnicity` (for each PC)
- **Outcome equation**: `birth_weight ~ c*GDM + Σ(b_i*PC_i) + sex_bw*sex + GA_bw*GA + ethnicity_bw*ethnicity`
- **Indirect effect**: `total_indirect := Σ(a_i × b_i)` across all PCs
- **Total effect**: `total := c + total_indirect` (direct + indirect)
- **Proportion mediated**: `prop_mediated := total_indirect / total`

**3. Sensitivity analysis**:
- Test PC numbers from 1 to 10 (GUSTO) or 1 to 15 (Gen3G)
- Select PCs explaining ≥75% cumulative variance
- GUSTO optimal PCs: 8 (GENCODEv45), 6 (GENCODE+), 5 (lr-assembly)
- Gen3G optimal PCs: 14 (GENCODEv45), 7 (GENCODE+), 6 (lr-assembly)

**4. Cross-tissue comparison**:
- GUSTO Placenta vs Adipose vs Fibroblasts assemblies
- Indirect effects and proportion mediated across tissue-specific annotations

**Gene Ontology enrichment analysis**:

**1. EnrichR databases**:
- GO_Biological_Process_2023
- GO_Molecular_Function_2023
- GO_Cellular_Component_2023
- KEGG_2021_Human
- WikiPathway_2021_Human
- Reactome_2022
- MSigDB_Hallmark_2020

**2. Gene universe definition**:
- **GENCODEv45 universe**: All genes in GENCODE v45 annotation
- **lr-assembly universe**: All ENSG-prefixed genes in filtered assembly
- Background-aware enrichment for accurate odds ratio calculation

**3. Ensembl-to-Symbol conversion**:
- org.Hs.eg.db for human gene annotation
- Remove version numbers from Ensembl IDs
- Filter out genes without valid HGNC symbols

**4. Enrichment workflows**:
- **Differential transcript expression (DTE)**: Map significant transcripts (padj < 0.1) to genes
- **Differential gene expression (DEG)**: Direct gene-level analysis
- **Top 50 isoform-rich genes**: Enrichment analysis on genes with highest transcript diversity
- **Mediating transcripts**: GO analysis on significant mediators (p < 0.05)

**5. Visualization**:
- Lollipop plots showing top 15 enriched terms
- Point size proportional to odds ratio
- Color by database source
- -log10(adjusted p-value) on x-axis
- Annotation labels for overlap and odds ratio

**Functional annotation of mediators**:

**1. GO Biological Process hierarchy**:
- Map genes to top-level GO BP terms
- GOBPPARENTS for hierarchical traversal
- Assign functional categories: Immune, Development, Regulation, Metabolism, Cell Organization, Unknown

**2. Category prioritization**:
- Priority order: Immune > Development > Regulation > Metabolism > Cell Organization > Unknown
- One category per gene (highest priority if multiple assignments)

**3. Gene symbol annotation**:
- HGNC symbols via org.Hs.eg.db
- Labeled visualization: "SYMBOL\n(transcript_id)"

**CSH1 transcript structure visualization**:

**1. Protein domain mapping**:
- HMMER Pfam domain predictions (E-value ≤ 1e-5)
- Protein coordinates mapped to genomic coordinates via SQANTI3 CDS annotations
- Domain merging for overlapping same-type domains

**2. ggtranscript plotting**:
- Exon-intron structure for all CSH1 isoforms
- Strand-aware intron arrows
- Protein domain tracks (SH, GHRL-hormone domains)
- Highlighting of significant mediators (bold transcript IDs, thicker lines)

**3. Structural categories**:
- Color by SQANTI3 classification (FSM, ISM, NIC)
- Y-axis ordering by structural category then transcript ID

**Output files**:

**Data files - Mediation analysis**:
- **all_results_tables.RData**: Comprehensive results (DE tables, mediation, VST data)
- **mediation_union_transcripts.RData**: Union set mediation results for cross-cohort validation
- **PCsensitivity.GUSTO.csv**: PC sensitivity indirect effects (GUSTO lr-assembly)
- **PCsensitivity.totals.GUSTO.csv**: PC sensitivity total effects (GUSTO lr-assembly)
- **PCsensitivity.GUSTO.gencode.csv**: PC sensitivity (GUSTO GENCODEv45)
- **PCsensitivity.totals.GUSTO.gencode.csv**: Total effects (GUSTO GENCODEv45)
- **PCsensitivity.GUSTO.gencode_plus.csv**: PC sensitivity (GUSTO GENCODE+)
- **PCsensitivity.totals.GUSTO.gencode_plus.csv**: Total effects (GUSTO GENCODE+)
- **PCsensitivity.GUSTO.adipose.csv**: PC sensitivity (GUSTO Adipose)
- **PCsensitivity.totals.GUSTO.adipose.csv**: Total effects (GUSTO Adipose)
- **PCsensitivity.GUSTO.fibroblasts.csv**: PC sensitivity (GUSTO Fibroblasts)
- **PCsensitivity.totals.GUSTO.fibroblasts.csv**: Total effects (GUSTO Fibroblasts)
- **PCsensitivity.Gen3G.csv** through **PCsensitivity.totals.Gen3G.gencode_plus.csv**: Gen3G PC sensitivity files

**Data files - Enrichment analysis**:
- **EnrichR_Results/EnrichR_Results_top50isorichgenes/results.csv**: Top 50 isoform-rich genes enrichment
- **top50isorichgenes_GOMF.csv**: GO Molecular Function for top 50 genes
- **EnrichR_Results/EnrichR_Results_DTE/results_GUSTO_DTE.csv**: GUSTO DTE enrichment summary
- **EnrichR_Results/EnrichR_Results_DTE/results_Gen3G_DTE.csv**: Gen3G DTE enrichment summary
- **EnrichR_Results/EnrichR_Results_DGE/results_GUSTO_DGE.csv**: GUSTO DEG enrichment summary
- **EnrichR_Results/EnrichR_Results_DGE/results_Gen3G_DGE.csv**: Gen3G DEG enrichment summary
- **EnrichR_Results/EnrichR_Results_DTE/\*.csv**: Individual database results for DTE
- **EnrichR_Results/EnrichR_Results_DGE/\*.csv**: Individual database results for DEG

**Figures - Mediation analysis**:
- **Figure 3F**: Direct and indirect effects across annotations (GUSTO only)
  - Total effect, indirect effect, and proportion mediated
  - Error bars showing 95% CI
  - Selected PCs explaining ≥75% variance
- **Figure 4A**: Mediation effects for union of top 10 significant mediators
  - Forest plot with error bars (95% CI)
  - Faceted by cohort (GUSTO, Gen3G)
  - Color by functional category
  - Bold gene symbols, filled points for p < 0.05
- **Figure 4B**: CSH1 isoform structure with protein domains
  - ggtranscript exon-intron plot
  - Pfam domain tracks (SH, GHRL-hormone)
  - Highlighted significant mediators (ESPRESSO:chr17:12985:70, ESPRESSO:chr17:12985:62)

**Figures - PC sensitivity analysis**:
- **Figure S8A**: PC mediation sensitivity for GUSTO (GENCODEv45, GENCODE+, lr-assembly)
  - Two-panel plot: Indirect Effect and % Variance Explained
  - X-axis: Number of PCs (1-10)
  - Horizontal reference lines at 0 (indirect) and 75% (variance)
- **Figure S8B**: PC mediation sensitivity for Gen3G (GENCODEv45, GENCODE+, lr-assembly)
  - Same layout as S8A but for Gen3G cohort
  - X-axis: Number of PCs (1-15)
- **Figure S8D**: PC mediation sensitivity for GUSTO tissues (Placenta, Adipose, Fibroblasts)
  - Cross-tissue comparison of mediation effects
  - Color scheme: seagreen (Placenta), #ff68a1 (Adipose), #e68613 (Fibroblasts)

**Figures - Enrichment analysis**:
- **EnrichR_Results/EnrichR_Results_DTE/Plots/\*_lollipop_DTE.pdf**: 
  - Lollipop plots for top 15 enriched terms (DTE)
  - Datasets: GUSTO_assembly_sig, GUSTO_gencode_sig, Gen3G_assembly_sig, Gen3G_gencode_sig
  - Point size = odds ratio, color = database source
- **EnrichR_Results/EnrichR_Results_DGE/Plots/\*_lollipop_DGE.pdf**:
  - Lollipop plots for top 15 enriched terms (DEG)
  - Same datasets as DTE analysis

To generate all mediation analysis results, PC sensitivity plots, enrichment analyses, and CSH1 structural visualization:
📜 [IVD_mediation_and_GOEA.R](IVD_mediation_and_GOEA.R)

---