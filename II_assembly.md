---
output:
  html_document: default
  pdf_document: default
---

# Long-read transcriptome assembly reveals vast isoform diversity in the placenta associated with metabolic and endocrine function

## II. Transcriptome assembly

**Overview**: This section describes the long-read transcriptome assembly pipeline using Oxford Nanopore Technology (ONT) sequencing data. The workflow processes raw nanopore signals through basecalling, orientation/trimming, error correction, alignment, and de novo transcript assembly to generate a comprehensive placental transcriptome.

**Key pipeline stages**:

1. **Basecalling**: Convert raw electrical signals to nucleotide sequences
2. **Orientation & trimming**: Identify full-length cDNA and remove adapters
3. **Error correction**: Use short-read data to correct sequencing errors
4. **Alignment**: Map corrected reads to reference genome
5. **Assembly**: Construct transcript models using ESPRESSO

### Requirements

- samtools v1.9
- pychopper v2.7.9
- fmlrc2 
- minimap2 v2.24
- ESPRESSO

### Setup

```bash
################################################################################
# Directory Structure and File Paths
################################################################################

# Primary data directories
DIR_RAW=/path/to/raw_reads                    # Raw ONT FASTQ files
DIR_TRIM=/path/to/trimmed_reads               # Pychopper-processed reads
DIR_ALIGN=/path/to/alignment_files            # Minimap2 alignments
DIR_INDEX=/path/to/preprocessing_indices_directory     # Reference indices
DIR_FMLRC2=/path/to/corrected_reads_directory          # Error-corrected reads
DIR_REPORTS=/path/to/reports_directory                 # QC reports
DIR_METRICS=/path/to/metrics_directory                 # Processing statistics
DIR_ESPRESSO=/path/to/txome_assembly                   # ESPRESSO assembly output
DIR_ESPRESSO_SCRIPTS=/path/to/ESPRESSO_git_repo_scripts # ESPRESSO utility scripts

# Reference files
GENOME_FASTA=/path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna  # hg38 genome
REFERENCE_GTF=/path/to/gencode_v45_annotation.gtf                      # GENCODE v45
SAMPLES=/path/to/samples.tsv_for_ESPRESSO                             # Sample metadata

# Sample identifiers
LIBIDs=() # SampleIDs as in .fastq file names, e.g. Sample01.fastq
```

### ONT basecalling

**Purpose**: Convert raw electrical signals (FAST5) from nanopore sequencing into nucleotide sequences (FASTQ). Basecalling accuracy directly impacts downstream assembly quality.

**Process**: FAST5 files of the raw ONT signals for each sample were processed with Guppy vX.X.X using the XXX.cfg configuration (Oxford Nanopore Technologies) for base calling on a GPU. FASTQ files labeled "pass" were concatenated together for each sample.

**Note**: The specific Guppy version and configuration should be documented for reproducibility, as basecalling models significantly impact accuracy.

### Long-read orientation & adapter trimming

**Purpose**: 

- **Identify full-length cDNA**: Distinguish complete transcripts from partial/degraded reads
- **Orient reads**: Ensure consistent 5' to 3' directionality 
- **Remove adapters**: Clean sequences from library preparation artifacts

**Tool**: Pychopper2 specifically designed for ONT cDNA library processing with PCS111 kit.

For each \$LIBID *n* in \$\{LIBIDs[1...*n*]\}

```bash
FASTQ=${DIR_RAW}/${LIBID}.fastq.gz

# Pychopper processing with comprehensive output tracking
pychopper -r ${DIR_REPORTS}/${LIBID}_pychopper_report.pdf \
  -u ${DIR_TRIM}/${LIBID}.unclassified.fastq \
  -w ${DIR_TRIM}/${LIBID}.rescued.fastq \
  -k PCS111 \
  -t ${THREADS} -S ${DIR_METRICS}/${LIBID}_pychopper_statistics.tsv \
  ${FASTQ} ${DIR_TRIM}/${LIBID}.fullLength.fastq
```

**Output categories**:

- **Full-length**: Complete cDNA with both 5' and 3' adapters detected
- **Rescued**: Reads with partial adapter detection but likely full-length
- **Unclassified**: Reads without clear adapter patterns (excluded from assembly)

### Long-read error correction

**Purpose**: Correct sequencing errors in ONT reads using high-accuracy short-read data. ONT reads have ~5-15% error rate, while short reads have ~0.1% error rate.

**Strategy**: FMLRC2 uses k-mer frequencies from short-read data to identify and correct likely sequencing errors while preserving genuine splice junctions and biological variation.

**Short-read source**: GUSTO villous placenta samples (N=200) provide comprehensive k-mer coverage for error correction.

For each \$LIBID *n* in \$\{LIBIDs[1...*n*]\}
```bash
# Combine full-length and rescued reads for maximum transcript recovery
cat ${DIR_TRIM}/${LIBID}.fullLength.fastq ${DIR_TRIM}/${LIBID}.rescued.fastq \
  > ${DIR_TRIM}/${LIBID}.fullLength.and.rescued.fastq
  
# Error correction using pre-built BWT index from short-read data
fmlrc2 -t ${THREADS} ${DIR_INDEX}/msbwt_GUSTO.npy \
  ${DIR_TRIM}/${LIBID}.fullLength.and.rescued.fastq \
  ${DIR_FMLRC2}/${LIBID}.fullLength.and.rescued.corrected.fasta
```

**Benefits of error correction**:

- Improved alignment accuracy and splice junction detection
- Reduced false positive isoforms from sequencing artifacts
- Better consensus building in transcript assembly

### Long-read alignment and filtering

**Purpose**: Map error-corrected long reads to the reference genome to identify exon-intron structures and splice patterns.

**Alignment strategy**: 

- **Minimap2**: Optimized for long-read spliced alignment
- **Splice-aware mode**: Detects splice junctions without strict GT-AG assumptions
- **Relaxed flanking**: `--splice-flank=no` accommodates non-canonical splice sites

**Quality control**: Remove secondary alignments and focus on primary chromosomes to reduce computational complexity and false discoveries.

For each \$LIBID *n* in \$\{LIBIDs[1...*n*]\}
```bash
# Minimap2 alignment with splice-aware parameters
# k=14, w=4: Optimized for transcript-length reads
# --secondary=no: Keep only best alignment per read
minimap2 -a -x splice --splice-flank=no \
         -w 4 -k 14 -t ${THREADS} -2 -K 2G --secondary=no \
         ${GENOME_FASTA} ${DIR_FMLRC2}/${LIBID}.fullLength.and.rescued.corrected.fasta \
         > ${DIR_ALIGN}/${LIBID}.sam
  
# Filter secondary alignments and coordinate sort
# FLAG 256 = secondary alignment
samtools view -@ ${THREADS} -F 256 -u ${DIR_ALIGN}/${LIBID}.sam \
   | samtools sort -@ ${THREADS} - -O BAM -o ${DIR_ALIGN}/${LIBID}.bam
samtools index ${DIR_ALIGN}/${LIBID}.bam

rm ${DIR_ALIGN}/${LIBID}.sam

# Subset to standard chromosomes only
# Excludes scaffolds, patches, and alternative haplotypes
chrs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
      chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM)

# Split alignments by chromosome for parallel processing
mkdir ${DIR_ALIGN}/split_by_chr/${LIBID}

for chr in "${chrs[@]}"
do
   samtools view -@ ${THREADS} -b ${DIR_ALIGN}/${LIBID}.bam ${chr} \
                  > ${DIR_ALIGN}/split_by_chr/${LIBID}_${chr}.bam
   samtools index ${DIR_ALIGN}/split_by_chr/${LIBID}/${LIBID}_${chr}.bam
done
```

**Rationale for chromosome splitting**:

- **Computational efficiency**: Parallel processing across chromosomes
- **Memory management**: Reduced RAM requirements for large datasets
- **Quality control**: Focus on well-annotated genomic regions

### Transcriptome assembly with ESPRESSO

**Purpose**: Construct de novo transcript models from aligned long reads. ESPRESSO provides optimal balance between transcript discovery and false positive control.

**ESPRESSO advantages**:

- **Error-aware**: Accounts for remaining sequencing errors
- **Reference-guided**: Uses existing annotations to improve accuracy
- **Isoform-centric**: Specifically designed for transcript assembly
- **Scalable**: Efficient processing of large datasets

**Assembly strategy**: Process each chromosome separately to manage computational complexity and enable parallel execution.

#### ESPRESSO_S (Sample processing)

**Purpose**: Process individual samples and prepare data for transcript clustering.

```bash
# ESPRESSO_S: Initial sample processing and read clustering
ESPRESSO=/path/to/scripts

perl ESPRESSO_S.pl -L ${SAMPLES} -F ${GENOME_FASTA} -A ${REFERENCE_GTF} \
                                 -O ${DIR_ESPRESSO} -T ${THREADS}
  
# Split large datasets for efficient processing
# Divides reads into manageable chunks for ESPRESSO_C
python3 ${DIR_ESPRESSO_SCRIPTS}/split_espresso_s_output_for_c.py \
   --orig-work-dir ${DIR_ESPRESSO} \
   --new-base-dir ${DIR_ESPRESSO}/split_S_for_C \
   --target-reads-per-c 2000000 \
   --genome-fasta ${GENOME_FASTA}
```

**Key functions**:

- **Read grouping**: Cluster reads by splice patterns
- **Initial filtering**: Remove low-quality alignments
- **Data partitioning**: Prepare for parallel processing

#### ESPRESSO_C (Correction/Consensus)

**Purpose**: Build consensus transcript models from clustered reads and correct remaining errors.

```bash
# Get updated sample list after splitting
SAMPLES_UPDATED=( $(awk '{print $1}' ${DIR_ESPRESSO}/samples.tsv.updated) )
```

For each \$SAMPLE *n* in \$\{SAMPLES_UPDATED[1...*n*]\}
```bash
# ESPRESSO_C: Consensus building and error correction
# -X 0: No additional filtering (already done in S step)
perl ESPRESSO_C.pl -I ${DIR_ESPRESSO}/split_S_for_C/${SAMPLE} -X 0 -T ${THREADS} \
                   -F ${DIR_ESPRESSO}/split_S_for_C/fastas/${SAMPLE}.fa 
```

**Process details**:

- **Consensus calling**: Build representative sequences from read clusters
- **Error correction**: Identify and fix remaining sequencing artifacts  
- **Splice junction refinement**: Improve splice site precision
- **Coverage assessment**: Calculate support levels for each transcript

#### ESPRESSO_Q (Quantification and finalization)

**Purpose**: Quantify transcript expression, merge samples, and generate final assembly.

```bash
cd ${DIR_ESPRESSO}

# Combine results from parallel C processing
python3 ${DIR_ESPRESSO_SCRIPTS}/combine_espresso_c_output_for_q.py \
                   --c-base-work-dir ${DIR_ESPRESSO}/split_S_for_C \
                   --new-base-dir ${DIR_ESPRESSO}/combine_C_for_Q
  
# ESPRESSO_Q: Final quantification and assembly generation
perl $ESPRESSO/ESPRESSO_Q.pl \
     -L ${DIR_ESPRESSO}/combine_C_for_Q/samples.tsv.updated \
     -A ${GTF} -T ${THREADS}
```

**Final outputs**:

- **GTF annotation**: Transcript coordinates and structure
- **FASTA sequences**: Transcript nucleotide sequences  
- **Expression matrix**: Read counts per transcript per sample
- **Quality metrics**: Assembly statistics and validation data

---

## Assembly Quality Considerations

### **Error Sources and Mitigation**

1. **Sequencing errors**: Addressed by FMLRC2 correction
2. **Alignment artifacts**: Minimized by splice-aware alignment
3. **Assembly errors**: Controlled by ESPRESSO's error-aware algorithms
4. **Reference bias**: Balanced by de novo discovery capabilities

### **Quality Control Metrics**

- **Read utilization**: Percentage of reads contributing to final assembly
- **Splice junction validation**: Consistency with short-read data
- **Expression correlation**: Agreement between samples
- **Reference concordance**: Overlap with known annotations

### **Expected Outcomes**

- **Novel isoforms**: Tissue-specific splice variants
- **Novel genes**: Previously unannotated transcriptional units
- **Improved annotations**: Extended UTRs and alternative promoters
- **Quantitative data**: Expression levels for all assembled transcripts

The final assembly provides a comprehensive catalog of placental transcript diversity, serving as the foundation for downstream functional and comparative analyses.

---