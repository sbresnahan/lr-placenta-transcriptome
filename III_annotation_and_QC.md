---
output:
  html_document: default
  pdf_document: default
---

# Long-read assembly reveals vast transcriptional complexity in the placenta associated with metabolic and endocrine function

## III. Transcriptome assembly annotation, QC and analysis

**Overview**: This section describes the annotation and quality control of the assembled placental transcriptome using SQANTI3. The workflow classifies transcript structures against reference annotations, validates transcript boundaries using orthogonal datasets, and applies comprehensive filtering to generate a high-confidence transcript catalog.

*Note: assembly QC analyses performed here were repeated for the GTEx v9 long-read data,
which are available under controlled access via dbGaP [accession phs000424.v9](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000424.v9) 
and on [AnVIL](https://anvil.terra.bio/#workspaces/anvil-datastorage/AnVIL_GTEx_V9_hg38).
See the manuscript and [orthogonal_for_GTEx.txt](Datasets/orthogonal_for_GTEx.txt) for details.*

**Objectives**:

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

**Purpose**: Generate comprehensive visualizations and statistical analyses of the long-read transcriptome assembly, comparing it to GENCODE v45 reference and GTEx v9 tissue atlases to characterize transcript diversity, structural categories, and transcriptional complexity.

**Analyses performed**:
- **RNA quality control**: RIN score comparison between GDM and control groups for ONT long-read and GUSTO short-read cohorts
- **Assembly statistics**: Gene and isoform counts compared to GENCODE v45
- **Structural classification**: Distribution of transcript categories before and after quality filtering
- **Transcript characteristics**: Length and exon count distributions by structural category
- **Multi-tissue comparison**: Transcriptional breadth (genes and isoforms detected) across GTEx v9 tissues and placenta
- **Expression saturation**: Isoform detection at varying expression thresholds across tissues
- **Isoform complexity**: Transcripts per gene analysis (assembly vs reference and across tissues)

**Output figures**:
- **RIN.png**: RIN score distributions for ONT and GUSTO cohorts
- **Figure 1B**: Gene and isoform counts (GENCODE v45 vs lr-assembly)
- **Figure 1C**: Structural category distribution (before vs after filtering)
- **Figure 1D**: Transcriptional breadth across GTEx tissues and placenta
- **Figure 1E**: Saturation curves showing isoforms detected at expression thresholds
- **Figure 1F**: Transcriptional complexity (isoforms per gene) across tissues
- **Figure S1A**: Transcript length distribution by structural category
- **Figure S1B**: Exon count distribution by structural category  
- **Figure S1C**: Forest plot of mean isoforms per gene (lr-assembly vs GENCODE v45)

To generate RIN quality control and Figures 1B-F and S1A-C:
📜 [IIIB_raw_assembly_analysis.R](IIIB_raw_assembly_analysis.R)

### C. LR-assembly Quality Control, Filtering, and Characterization

**Purpose**: Apply multi-layered quality control filters using orthogonal evidence from short-read data, CAGE-seq, poly(A) sites, and artifact detection to remove low-confidence transcripts and technical artifacts, then comprehensively characterize the high-quality filtered assembly.

**Quality control strategy**:

**1. Orthogonal evidence filtering**:
- **5' end support (TSS)**: CAGE-seq peaks (±5kb) or reference TSS proximity (±100bp)
- **3' end support (TTS)**: Poly(A) site annotations (±100bp), poly(A) motifs (-41 to 0bp), or reference TTS proximity (±100bp)
- **Splice junction support**: Short-read validation from GUSTO cohort (n=200) with minimum 3 reads per junction
- **Full-length support**: Long-read evidence from ONT data (n=72) with minimum 3 full-length reads
- **TSS ratio calculation**: First exon reads / last exon reads from short-read data to detect 5' degradation

**2. Artifact detection and removal**:
- **Poly(A) internal priming**: >80% A-content downstream of TTS (excludes FSM transcripts)
- **RT-switching artifacts**: Template switching during reverse transcription (applied to NNC, Genic, Antisense, Fusion, Intergenic categories only)
- **Degradation products**: TSS ratio <0.5 or predicted NMD targets (applied to novel categories only)

**3. Post-filtering characterization**:
- **Structural analysis**: Length, exon count, and complexity distributions after quality filtering
- **Novelty classification**: Known vs novel at both gene and isoform levels
- **Tissue enrichment**: TissueEnrich analysis on top 100-4000 expressed genes with odds ratio calculations
- **Isoform complexity**: Comparison of transcripts per gene between filtered assembly and GENCODE v45
- **Protein domain mapping**: Pfam domain visualization using HMMER predictions
- **Coding sequence validation**: 
  - BLASTP against placenta-specific MS/MS database
  - BLASTP against full UniProt database
  - CPC2 coding probability scores
  - TransDecoder ORF predictions
  - Unique peptide identification per isoform

**4. Transcript feature analysis**:
- Expression levels (FL reads) by novelty category
- Length and exon count distributions by novelty
- Correlation between gene complexity and exon counts
- GO/pathway enrichment of high-complexity genes

**Output components**:

**Filtered assembly files**:
- **classification_filtered.txt**: SQANTI3 classification for high-quality transcripts
- **corrected_filtered.gtf**: GTF annotation with filtered transcripts
- **corrected_filtered.fasta**: Nucleotide sequences
- **corrected_filtered.faa**: Protein sequences
- **junctions_filtered.txt**: Splice junction annotations
- **filter_reasons.txt**: Documentation of why each transcript was removed
- **tx2g.RData**: Transcript-to-gene mapping objects

**Quality control figures (Pre-filtering metrics)**:
- **Figure S2A**: Full-length read distribution by structural category (boxplot with filter threshold)
- **Figure S2B**: Short-read splice junction support minimum coverage (boxplot with filter threshold)
- **Figure S2C**: Distance from TSS to CAGE/ATAC peak (histogram with ±5kb thresholds)
- **Figure S2D**: Distance from TTS to poly(A) site (histogram with ±100bp thresholds)
- **Figure S2E**: Distance from TTS to poly(A) motif (histogram with -41 to 0bp thresholds)
- **Figure S2F**: Sequencing artifact detection by category (faceted bar plot: poly(A) intrapriming, RT switching, NMD/degradation)

**Filtered assembly characterization**:
- **Figure S3A**: Transcript length distribution post-filtering (boxplot by structural category)
- **Figure S3B**: Exon count distribution post-filtering (stacked bar chart)
- **Figure S3C**: Known vs novel gene and isoform classification (boxplot)
- **Figure S3D**: Tissue-specific enrichment across top 100-2000 expressed genes (jitter plot with odds ratios, placenta highlighted in black)

**Transcript features**:
- **Figure S4A**: Transcript abundance (FL reads) by novelty category (within known genes; known vs novel genes)
- **Figure S4B**: Transcript length by novelty category (boxplot comparisons)
- **Figure S4C**: Exon count by novelty category (boxplot comparisons)
- **Figure S4D**: Correlation between transcripts per gene and total exons (all transcripts)
- **Figure S4E**: Correlation between transcripts per gene and total exons (highly-expressed: log10 TPM > 2.5)

**Coding potential**:
- **Figure S5A**: Percentage of isoforms with ORFs by structural category
- **Figure S5B**: ORF length distribution by structural category (TransDecoder predictions)
- **Figure S5C**: Coding probability distribution by structural category (CPC2 scores)

**Main manuscript figures**:
- **Figure 2A**: Isoform complexity scatter plot (filtered assembly vs GENCODE v45) with marginal histograms, highlighting CSH1 and GH2 genes
- **Figure 2B**: CSH1 gene structure visualization with Pfam protein domains (SH/GHRL-hormone domain) using ggtranscript
- **Figure 2C**: CDS support validation - percentage of coding isoforms (>50% coding probability) with unique peptide support >95% identity from placenta MS/MS (solid bars) vs full UniProt (transparent bars)

To generate filtered assembly files and all quality control, characterization, and validation figures:
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