---
output:
  html_document: default
  pdf_document: default
---

# Long-read transcriptome assembly reveals vast isoform diversity in the placenta associated with metabolic and endocrine function

## I. Preparation of datasets for assembly

**Overview**: This section prepares the orthogonal datasets required for long-read transcriptome assembly validation and quality control. The pipeline integrates multiple data types to provide comprehensive support for transcript structure validation:

- **Short-read RNA-seq** (GUSTO, N=200): Provides splice junction validation and read coverage support
- **CAGE-seq** (FANTOM5): Identifies authentic transcription start sites (TSS)
- **DNase-seq** (ENCODE): Maps open chromatin regions associated with active promoters
- **Poly(A) sites** (APASdb + PolyA_DB): Validates transcription termination sites (TTS)

Short-read RNA-seq from term, villous placenta from GUSTO (N = 200) were 
aligned to hg38 using STAR with settings recommended by ENCODE. Splice junctions 
with < 10 supporting reads were filtered. CAGE-seq reads from placenta tissue 
(N = 1) and cells (pericyte, N = 2; epithelial, N = 4, trophoblast, N = 1) were 
retrieved from the FANTOM5 Project, trimmed with cutadapt v4.1 to remove linker 
adapters, EcoP15 and 5'poly-G sequences, and filtered reads were aligned to h38 
using STAR with settings recommended by ENCODE. CAGE tags were counted and 
clustered using CAGEr in R. DNAse-seq peaks mapped to h38 were retrieved in BED 
format from the ENCODE consortium and intersected with bedtools. SAPAS peaks in 
placenta tissue mapped to hg38 were retrieved in BED format from APASdb, and 
human polyA motifs annotated in PolyA_DB were retrieved from PolyA-miner.

### Setup

```bash
################################################################################
# Directory Structure and File Paths
################################################################################

# Primary data directories
DIR_RAW=/path/to/raw_reads                    # Raw sequencing data
DIR_TRIM=/path/to/trimmed_reads               # Quality-trimmed reads
DIR_ALIGN=/path/to/alignment_files            # STAR alignment outputs
DIR_CAGE=/path/to/processed_cage-seq_data     # CAGE-seq processing results
DIR_DNAse=/path/to/ENCODE_DNase-seq_narrowPeaks  # ENCODE DNase-seq peaks
DIR_ORTHOGONAL=/path/to/orthogonal_datasets_for_sqanti3  # Final SQANTI3 inputs

# Processing and temporary directories
DIR_INDEX=/path/to/preprocessing_indices_directory    # Reference genome indices
DIR_FMLRC2=/path/to/corrected_reads_directory         # FMLRC2 error correction
DIR_REPORTS=/path/to/reports_directory                # QC and processing reports
DIR_METRICS=/path/to/metrics_directory                # Alignment and trimming metrics

# Reference genome and annotation files
GENOME_FASTA=/path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna  # hg38 genome
REFERENCE_GTF=/path/to/gencode_v45_annotation.gtf                      # GENCODE v45
```

### A. Prepare reference indices for STAR alignment

**Purpose**: Generate STAR genome indices optimized for RNA-seq alignment to hg38 reference genome.

#### Requirements

- star v2.7.4a

#### STAR genome index for hg38

```bash
# Generate STAR genome index for hg38 reference
# Uses default parameters optimized for human genome size and complexity
STAR --runThreadN ${THREADS} --runMode genomeGenerate \
     --genomeDir ${DIR_INDEX}/star-2.7.4a_GCA_000001405.15_GRCh38_no_alt_analysis_set \
     --genomeFastaFiles ${GENOME_FASTA}
     
GENOME_INDEX_STAR=${DIR_INDEX}/star-2.7.4a_GCA_000001405.15_GRCh38_no_alt_analysis_set
```

### B. GUSTO short-read alignments to reference genome

**Purpose**: Process short-read RNA-seq data from 200 placental samples to generate:

1. High-quality splice junction annotations for long-read validation
2. Read coverage data for transcript support validation
3. Error-corrected reads for FMLRC2 long-read polishing

#### Requirements

- samtools v1.9
- fastp v0.23.2
- star v2.7.4a

#### Read trimming & quality control

```bash
# Define sample identifiers
LIBIDs=() # SampleIDs as in .fastq file names, e.g. Sample01.fastq.gz
```

For each \$LIBID *n* in \$\{LIBIDs[1...*n*]\}
```bash
# Trim adapters and low-quality bases using fastp
# Automatic adapter detection and quality filtering
fastp \
   -i ${DIR_RAW}/${LIBID}_1.fastq.gz -o ${DIR_TRIM}/${LIBID}_1.fastq.gz \
   -I ${DIR_RAW}/${LIBID}_2.fastq.gz -O ${DIR_TRIM}/${LIBID}_2.fastq.gz \
   -w ${THREADS} 
```

#### Read alignment with STAR in two-pass mode

**Two-pass mode benefits**:
- First pass: Discovers novel splice junctions
- Second pass: Re-aligns reads using discovered junctions for improved accuracy

For each \$LIBID *n* in \$\{LIBIDs[1...*n*]\}
```bash
# STAR alignment with ENCODE-recommended parameters
# Optimized for splice junction discovery and accurate alignment
STAR --genomeDir ${GENOME_INDEX_STAR} \
     --readFilesIn ${DIR_TRIM}/${LIBID}_1.fastq.gz ${DIR_TRIM}/${LIBID}_2.fastq.gz \
     --readFilesCommand zcat --runThreadN ${THREADS} \
     --genomeLoad NoSharedMemory --outFilterMultimapNmax 20 \
     --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
     --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 \
     --alignIntronMin 20 --alignIntronMax 1000000 \
     --alignMatesGapMax 1000000 --outSAMheaderHD @HD VN:1.4 SO:coordinate \
     --outSAMunmapped Within --outFilterType BySJout \
     --outSAMattributes NH HI AS NM MD --outSAMtype BAM SortedByCoordinate \
     --sjdbScore 1 --outTmpDir ${DIR_TEMP}/${LIBID} \
     --outFileNamePrefix ${DIR_ALIGN}/${LIBID} \
     --outBAMsortingBinsN 200 \
     --limitBAMsortRAM 80000000000 \
     --twopassMode Basic
     
# Index BAM files for efficient downstream processing
samtools index ${DIR_ALIGN}/${LIBID}Aligned.sortedByCoord.out.bam
```

Extract high-quality aligned reads for FMLRC2 error correction:

```bash
# Extract only properly mapped reads (FLAG 4 = unmapped)
# Convert to FASTQ for FMLRC2 BWT construction
samtools fastq -F 4 -@ 12 ${DIR_ALIGN}/GUSTO/${LIBID}Aligned.sortedByCoord.out.bam \
               -1 ${DIR_ALIGN}/GUSTO/${LIBID}_1.fastq.gz \
               -2 ${DIR_ALIGN}/GUSTO/${LIBID}_2.fastq.gz
```

#### Prepare splice junction coverage file

**Purpose**: Aggregate splice junction support across all samples for transcript validation.

📜 Performed in R with [merge_SJs.R](merge_SJs.R).

### C. FMLRC2 - Long-read error correction preparation

**Purpose**: Build Burrows-Wheeler Transform (BWT) index from short-read data for long-read error correction. FMLRC2 uses short-read k-mers to correct sequencing errors in long reads while preserving splice junction information.

#### Requirements

- samtools v1.9
- [fmlrc2 v0.1.8](https://github.com/HudsonAlpha/fmlrc2/tree/v0.1.8)
- [ropebwt2](https://github.com/lh3/ropebwt2)

#### Build BWT index from short-read data

```bash
# Extract all short-read sequences for BWT construction
echo "Decompressing, combining, and isolating short reads"
gunzip -c ${DIR_ALIGN}/GUSTO/*.fastq.gz | \
   awk 'NR % 4 == 2' | \
   tr NT TN > ${DIR_FMLRC2}/TEMP/short-reads.seq

# Build compressed BWT index for efficient k-mer lookup
echo "Building BWT"
ropebwt3 build -t ${THREADS} -L -R ${DIR_FMLRC2}/TEMP/short-reads.seq | \
  tr NT TN | \
  fmlrc2-convert ${DIR_INDEX}/msbwt_GUSTO.npy
```

### D. Prepare CAGE-seq peaks for TSS validation

**Purpose**: Process CAGE-seq data to identify authentic transcription start sites. CAGE (Cap Analysis Gene Expression) specifically captures the 5' ends of capped mRNAs, providing precise TSS annotations.

**Data sources**: FANTOM5 Project placenta samples
- Tissue: N=1 
- Cell types: pericyte (N=2), epithelial (N=4), trophoblast (N=1)

**Technical note**: DRR accessions were uploaded without quality scores, requiring separate processing from SRR accessions.

#### Requirements

- cutadapt v4.1
- samtools v1.9
- star v2.7.4a
- bedtools v2.31.1
- sra-tools

#### Retrieve CAGE-seq data from FANTOM5

```bash
# FANTOM5 placenta CAGE-seq accessions
CAGE_IDS=(DRR009582 SRR530926 SRR530927 SRR530924 SRR530925 \
          DRR009276 DRR009277 DRR009278 DRR008721)
  
for ACCESSION in "${CAGE_IDS[@]}"
do
  prefetch $ACCESSION
  fasterq-dump $ACCESSION
done
```

#### Process DRR accessions (no quality scores)

```bash
DRR_IDS=(DRR009582 DRR009276 DRR009277 DRR009278 DRR008721)
```

For each \$FILE *n* in \$\{DRR_IDs[1...*n*]\}
```bash
# Convert FASTQ to FASTA (no quality scores available)
sed -n '1~4s/^@/>/p;2~4p' ${DIR_RAW}/${FILE}.fastq > ${DIR_RAW}/${FILE}.fasta

# Trim EcoP15 linker sequences (CAGE-seq library preparation artifact)
# Pattern: CAGCAG...TCGTATGCCGTCTTC (variable length middle section)
cutadapt -a ^CAGCAG...TCGTATGCCGTCTTC \
         --match-read-wildcards --minimum-length 15 \
         --cores=${THREADS} \
         -o ${DIR_TRIM}/${FILE}_atrimmed.fasta ${DIR_RAW}/${FILE}.fasta \
          > ${DIR_METRICS}/${FILE}_atrimmed.metrics.txt
  
rm ${DIR_RAW}/${FILE}.fasta

# Trim 5'poly-G sequences (common CAGE-seq artifact)
# These can interfere with accurate TSS mapping
cutadapt -g ^G -e 0 --match-read-wildcards \
         --cores=${THREADS} \
         -o TRIM/${FILE}_gtrimmed.fasta TRIM/${FILE}_atrimmed.fasta \
         > METRICS/${FILE}_gtrimmed.metrics.txt 
```

#### Process SRR accessions (with quality scores)

```bash
SRR_IDS=(SRR530926 SRR530927 SRR530924 SRR530925)
```

For each \$FILE *n* in \$\{SRR_IDs[1...*n*]\}
```bash
# Trim EcoP15 linker sequences
cutadapt -a ^CAGCAG...TCGTATGCCGTCTTC \
         --match-read-wildcards --minimum-length 15 \
         --cores=${THREADS} \
         -o ${DIR_TRIM}/${FILE}_atrimmed.fastq ${DIR_RAW}/${FILE}.fastq \
          > ${DIR_METRICS}/${FILE}_atrimmed.metrics.txt

# Trim 5'poly-G sequences
cutadapt -g ^G -e 0 --match-read-wildcards \
         --cores=${THREADS} \
         -o ${DIR_TRIM}/${FILE}_gtrimmed.fastq ${DIR_TRIM}/${FILE}_atrimmed.fastq \
          > ${DIR_METRICS}/${FILE}_gtrimmed.metrics.txt 
```

#### Align CAGE-seq reads to reference genome

For each \$FILE *n* in \$\{FILES[1...*n*]\}
```bash
# STAR alignment optimized for short CAGE tags
# More permissive settings due to short read length (15-50bp typical)
STAR --genomeDir ${GENOME_INDEX_STAR} \
     --readFilesIn ${DIR_TRIM}/${FILE}_gtrimmed.fasta \
     --runThreadN ${THREADS} \
     --outSAMtype BAM SortedByCoordinate \
     --genomeLoad NoSharedMemory --outFilterMultimapNmax 20 \
     --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
     --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 \
     --alignIntronMin 20 --alignIntronMax 1000000 \
     --alignMatesGapMax 1000000 --outSAMheaderHD @HD VN:1.4 SO:coordinate \
     --outSAMunmapped Within --outFilterType BySJout \
     --outSAMattributes NH HI AS NM MD \
     --sjdbScore 1 \
     --outBAMsortingBinsN 200 \
     --limitBAMsortRAM 80000000000 \
     --twopassMode Basic \
     --outTmpDir ${DIR_CAGE}/TEMP/${FILE} \
     --outFileNamePrefix ${DIR_CAGE}/${FILE}
```

#### Generate CAGE Tag Starting Sites (CTSS)

**Purpose**: Convert aligned CAGE reads to precise TSS coordinates with read count support.

For each \$FILE *n* in \$\{FILES[1...*n*]\}
```bash
# Filter for high-quality alignments only
# Q30 = 99.9% accuracy, essential for precise TSS mapping
samtools view -@ 8 -F 4 -q 30 -u ${DIR_CAGE}/${FILE}Aligned.sortedByCoord.out.bam \
  > ${DIR_CAGE}/${FILE}Aligned.sortedByCoord.out_filtered.bam

rm ${DIR_CAGE}/${FILE}Aligned.sortedByCoord.out.bam

# Convert BAM to BED format for easier processing
bamToBed -i ${DIR_CAGE}/${FILE}Aligned.sortedByCoord.out_filtered.bam \ 
  > ${DIR_CAGE}/${FILE}Aligned.sortedByCoord.out_filtered.bed

# Generate CTSS on positive strand
# Use 5' end of read (start position) as TSS
mkdir ${DIR_CAGE}/CTSS

awk 'BEGIN{OFS="\t"}{if($6=="+"){print $1,$2,$5}}' \
  ${DIR_CAGE}/${FILE}Aligned.sortedByCoord.out_filtered.bed \
  | sort -k1,1 -k2,2n \
  | groupBy -i stdin -g 1,2 -c 3 -o count \
  | awk -v x="$FILE" 'BEGIN{OFS="\t"}{print $1,$2,$2+1,  x  ,$3,"+"}' >> \
    ${DIR_CAGE}/CTSS/${FILE}.pos_ctss.bed

# Generate CTSS on negative strand
# Use 3' end of read (end position) as TSS for antisense transcripts
awk 'BEGIN{OFS="\t"}{if($6=="-"){print $1,$3,$5}}' \
  ${DIR_CAGE}/${FILE}Aligned.sortedByCoord.out_filtered.bed \
  | sort -k1,1 -k2,2n \
  | groupBy -i stdin -g 1,2 -c 3 -o count \
  | awk -v x="$FILE" 'BEGIN{OFS="\t"}{print $1,$2-1,$2, x  ,$3,"-"}' >> \
    ${DIR_CAGE}/CTSS/${FILE}.neg_ctss.bed

# Combine both strands and sort by genomic coordinates
cat ${DIR_CAGE}/CTSS/${FILE}.pos_ctss.bed ${DIR_CAGE}/CTSS/${FILE}.neg_ctss.bed \
  | sort -k1,1 -k2,2n > ${DIR_CAGE}/CTSS/${FILE}.ctss.bed

# Clean up intermediate files
rm ${DIR_CAGE}/${FILE}Aligned.sortedByCoord.out_filtered.bed
${DIR_CAGE}/CTSS/${FILE}.pos_ctss.bed
${DIR_CAGE}/CTSS/${FILE}.neg_ctss.bed
```

#### Aggregate CAGE tags across samples

**Purpose**: Combine CTSS from all samples and cluster nearby tags into coherent TSS regions.

Performed with CAGEr in R to generate `${DIR_CAGE}/FANTOM5_placenta_CAGE-seq.bed`.

### E. Prepare comprehensive TSS support (DNase-seq + CAGE-seq peaks)

**Purpose**: Integrate multiple lines of evidence for active transcription start sites:

- **CAGE-seq**: Direct measurement of capped mRNA 5' ends
- **DNase-seq**: Open chromatin regions indicating accessible promoters

**Rationale**: Combining both datasets increases sensitivity for TSS detection while maintaining specificity.

#### Data source
[DNase-seq peaks of placental tissue (N = 22) retrieved from ENCODE](https://www.encodeproject.org/search/?type=Experiment&control_type!=*&status=released&perturbed=false&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&biosample_ontology.organ_slims=placenta&biosample_ontology.term_name=placenta&assay_title=DNase-seq&biosample_ontology.classification=tissue&files.file_type=bed+narrowPeak)

```bash
# Merge DNase-seq peaks across all 22 placental samples
# Remove overlapping regions to create non-redundant peak set
cat ${DIR_DNAse}/*.bed \
  | bedtools sort -i - \
  | bedtools merge -i - > ${DIR_DNAse}/ENCODE_placenta_DNase-seq.bed
  
# Combine CAGE-seq & DNase-seq peaks for comprehensive TSS annotation
# This creates a union of both evidence types
cat ${DIR_CAGE}/FANTOM5_placenta_CAGE-seq.bed \
   ${DIR_DNAse}/ENCODE_placenta_DNase-seq.bed \
  | bedtools sort -i - \
  | bedtools merge -i - > ${DIR_ORTHOGONAL}/placenta_TSS.bed
```

### F. Prepare comprehensive TTS support (Poly(A) sites + motifs)

**Purpose**: Integrate multiple sources of transcription termination site evidence:

- **SAPAS**: Sequencing-based poly(A) site identification from tissue samples
- **PolyA_DB**: Curated database of poly(A) motifs and cleavage sites

**Rationale**: Combining experimental and curated data provides comprehensive TTS annotation.

#### Data sources
- [SAPAS peaks in placenta tissue were retrieved from APASdb](http://mosas.sysu.edu.cn/utr)
- [Human polyA motifs annotated in PolyA_DB were retrieved from PolyA-miner](https://github.com/YalamanchiliLab/PolyA-miner/blob/master/Human_hg38.PolyA_DB.bed)

```bash
# Combine SAPAS experimental data with curated poly(A) motifs
# Creates comprehensive TTS annotation for transcript 3' end validation
cat ${DIR_ORTHOGONAL}/APASdb_placenta_SAPAS.bed \
    ${DIR_ORTHOGONAL}/PolyA_DB_motifs.bed \
  | bedtools sort -i - \
  | bedtools merge -i - > ${DIR_ORTHOGONAL}/placenta_TTS.bed
```

---