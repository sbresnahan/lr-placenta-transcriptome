---
output:
  html_document: default
  pdf_document: default
---

# Long-read transcriptome assembly reveals vast isoform diversity in the placenta associated with metabolic and endocrine function

## III. Transcriptome assembly annotation, QC and analysis

**Overview**: This section describes the annotation and quality control of the assembled placental transcriptome using SQANTI3. The workflow classifies transcript structures against reference annotations, validates transcript boundaries using orthogonal datasets, and applies comprehensive filtering to generate a high-confidence transcript catalog.

**Key objectives**:

- **Structural classification**: Compare assembled transcripts to GENCODE v45 reference
- **Boundary validation**: Verify transcript start/end sites using CAGE-seq and poly(A) data
- **Quality assessment**: Identify and filter potential sequencing artifacts
- **Comprehensive characterization**: Generate detailed assembly statistics and visualizations

### Requirements

- [SQANTI3 v5.4](https://github.com/ConesaLab/SQANTI3/tree/v5.4)

### Setup

```bash
################################################################################
# Directory Structure and File Paths
################################################################################

# Assembly directories
DIR_ESPRESSO=/path/to/txome_assembly                    # ESPRESSO assembly output
DIR_ASSEMBLY=/path/to/txome_assembly_annotated          # SQANTI3 annotation output
DIR_SQANTI3_SCRIPTS=/path/to/SQANTI3_git_repo           # SQANTI3 installation
DIR_ORTHOGONAL=/path/to/orthogonal_datasets_for_sqanti3 # Validation datasets

# Reference files
GENOME_FASTA=/path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna  # hg38 genome
REFERENCE_GTF=/path/to/gencode_v45_annotation.gtf                      # GENCODE v45
```

### A. Isoform classification against GENCODEv45 via SQANTI3

**Purpose**: SQANTI3 performs comprehensive structural classification of assembled transcripts by comparing them to reference annotations. This process identifies the relationship between novel isoforms and known gene models, validates transcript boundaries using orthogonal datasets, and flags potential sequencing artifacts.

**Classification categories**:

- **Full Splice Match (FSM)**: Perfect match to reference transcript
- **Incomplete Splice Match (ISM)**: Subset of reference exons
- **Novel In Catalog (NIC)**: Novel combination of known splice sites
- **Novel Not in Catalog (NNC)**: Contains novel splice sites
- **Genic**: Overlaps gene but different structure
- **Antisense**: Antisense to annotated gene
- **Intergenic**: No overlap with annotated genes

**Validation datasets integration**:

- **CAGE peaks**: Validate transcription start sites
- **Poly(A) sites**: Validate transcription termination sites  
- **Short-read splice junctions**: Confirm splice site authenticity
- **Expression data**: Assess transcript support levels

```bash
################################################################################
# Prepare SQANTI3 Input Files
################################################################################

# Validation datasets prepared in previous sections
CAGE_PEAKS=${DIR_ORTHOGONAL}/placenta_TSS.bed                    # TSS validation
ASSEMBLY_GTF=${DIR_ESPRESSO}/combine_C_for_Q/sample_transcripts.gtf   # ESPRESSO assembly
LONG_READ_EXPRESSION=${DIR_ESPRESSO}/combine_C_for_Q/sample_quants.tsv # Expression matrix
SJOUT=${DIR_ORTHOGONAL}/GUSTO_filtered_SJ.out.tab                # Short-read splice junctions

################################################################################
# Run SQANTI3 Classification and Quality Control
################################################################################

# Comprehensive isoform classification with orthogonal validation
python3 ${DIR_SQANTI3_SCRIPTS}/sqanti3_qc.py \
  --force_id_ignore \
  -o PLACENTA_ASSEMBLY \
  -d ${DIR_ASSEMBLY} \
  --report html \
  --CAGE_peak ${CAGE_PEAKS} \
  --polyA_motif_list ${DIR_SQANTI3_SCRIPTS}/data/polyA_motifs/mouse_and_human.polyA_motif.txt \
  --polyA_peak ${DIR_ORTHOGONAL}/placenta_TTS.bed \
  --fl_count ${LONG_READ_EXPRESSION} \
  --ratio_TSS_metric mean \
  -c ${SJOUT} \
  ${ASSEMBLY_GTF} ${REFERENCE_GTF} ${GENOME_FASTA}
```

**Parameter explanations**:

- `--force_id_ignore`: Ignore transcript ID conflicts between assembly and reference
- `--report html`: Generate comprehensive HTML quality report
- `--CAGE_peak`: Integrate CAGE-seq data for TSS validation
- `--polyA_motif_list`: Use curated poly(A) motif database for TTS validation
- `--polyA_peak`: Integrate experimental poly(A) site data
- `--fl_count`: Include full-length read counts for expression-based filtering
- `--ratio_TSS_metric mean`: Use mean TSS ratio for degradation assessment
- `-c`: Include short-read splice junction support data

**SQANTI3 output files**:

- **Classification table**: Detailed structural annotations for each transcript
- **HTML report**: Interactive quality control dashboard
- **Junction analysis**: Splice site validation and novelty assessment
- **Expression integration**: Read support and tissue-specificity metrics
- **Artifact detection**: Identification of potential sequencing artifacts

### B. LR-assembly Analysis

**Purpose**: Generate comprehensive visualizations and statistical analyses of the raw (unfiltered) assembly to characterize transcript diversity, structural categories, and comparison with reference annotations.

**Analyses performed**:

- **Assembly statistics**: Gene and isoform counts compared to GENCODE v45
- **Structural distribution**: Proportion of each transcript category
- **Length and complexity**: Transcript length and exon count distributions
- **Multi-tissue comparison**: Assembly characteristics vs GTEx v9 tissues
- **Isoform complexity**: Transcripts per gene analysis

**Output figures**:

- **Figure 1B**: Gene and isoform counts (assembly vs reference)
- **Figure 1C**: Structural category distribution
- **Figure 1D**: Multi-tissue gene and isoform detection comparison
- **Figure 1E**: Transcriptional complexity across tissues
- **Figure S1A**: Transcript length distribution by structural category
- **Figure S1B**: Exon count distribution by structural category  
- **Figure S1C**: Isoform complexity comparison (assembly vs reference)

To generate Figures 1B-E and S1A-C:

📜 [IIIB_raw_assembly_analysis.R](IIIB_raw_assembly_analysis.R)

### C. LR-assembly Quality Control, Filtering, and Characterization

**Purpose**: Apply multi-layered quality control filters to remove potential artifacts and low-confidence transcripts, then characterize the filtered high-quality assembly.

**Quality control strategy**:

**1. Orthogonal evidence filtering**:

- **TSS support**: CAGE-seq or reference TSS proximity
- **TTS support**: Poly(A) sites or reference TTS proximity  
- **Splice junction support**: Short-read validation (≥3 reads)
- **Full-length support**: Long-read evidence (≥3 reads)

**2. Artifact detection and removal**:

- **Poly(A) internal priming**: Excessive A-content downstream of TTS
- **RT-switching artifacts**: Template switching during reverse transcription
- **Degradation products**: Low TSS ratio or predicted NMD targets

**3. Post-filtering characterization**:

- **Structural analysis**: Length and exon distributions after filtering
- **Novelty assessment**: Known vs novel gene and isoform classification
- **Tissue enrichment**: Functional relevance of top expressed transcripts
- **Protein validation**: Coding sequence support analysis

**Output components**:

**Filtered assembly files**:

- **GTF annotation**: High-confidence transcript coordinates
- **FASTA sequences**: Nucleotide and protein sequences
- **Expression matrix**: Quantification data
- **Filter documentation**: Reasons for transcript removal

**Quality control figures (Pre-filtering)**:

- **Figure S2A**: Full-length read distribution by structural category
- **Figure S2B**: Short-read splice junction coverage
- **Figure S2C**: CAGE/ATAC-seq TSS validation distances
- **Figure S2D**: Poly(A) site validation distances
- **Figure S2E**: Poly(A) motif distances  
- **Figure S2F**: Sequencing artifact detection summary

**Filtered assembly characterization**:

- **Figure S3A**: Post-filtering transcript length distribution
- **Figure S3B**: Post-filtering exon count distribution
- **Figure S3C**: Known vs novel gene/isoform classification
- **Figure S3D**: Tissue-specific enrichment of top expressed transcripts

**Main manuscript figures**:

- **Figure 2A**: Isoform complexity comparison (filtered assembly vs reference)
- **Figure 2B**: CSH1 gene structure visualization
- **Figure 2C**: Protein coding sequence support analysis

To generate filtered assembly files and Figures S2A-F, S3A-D, and 2A-C:

📜 [IIIC_assembly_QC_and_filtering.R](IIIC_assembly_QC_and_filtering.R)

---

## Quality Control Summary

### **Assembly Quality Metrics**

- **Classification accuracy**: Proportion of high-confidence transcript categories
- **Orthogonal validation**: Percentage with TSS/TTS/splice junction support
- **Artifact detection**: Identification and removal of technical artifacts
- **Expression validation**: Consistency across biological replicates

### **Filtering Impact**

- **Transcript retention**: Percentage of transcripts passing quality filters
- **Category enrichment**: Relative enrichment of high-confidence categories
- **Novel discovery**: Balance between sensitivity and specificity
- **Functional relevance**: Tissue-specific expression patterns

### **Final Assembly Characteristics**

- **Comprehensive coverage**: Representation of placental gene expression
- **High confidence**: Multiple lines of orthogonal evidence
- **Functional annotation**: Integration with protein databases
- **Comparative context**: Relationship to reference annotations and other tissues

The resulting filtered assembly provides a high-quality catalog of placental transcript diversity, suitable for downstream functional analyses, comparative studies, and biological interpretation.

---