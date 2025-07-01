# Long-read transcriptome assembly reveals vast isoform diversity in the placenta associated with metabolic and endocrine function

This repository contains the complete codebase and documentation for the analyses presented in our manuscript:

> **Long-read transcriptome assembly reveals vast isoform diversity in the placenta associated with metabolic and endocrine function**  
> Sean T Bresnahan, Hannah Yong, William H Wu, Sierra Lopez, Jerry Kok Yen Chan, Frédérique White, Pierre-Étienne Jacques, Marie-France Hivert, Shiao-Yng Chan, Michael I Love, Jonathan Y Huang, Arjun Bhattacharya  
> *bioRxiv*, 2025.06.26.661362  
> [https://doi.org/10.1101/2025.06.26.661362](https://doi.org/10.1101/2025.06.26.661362)

The pipeline integrates long-read sequencing, orthogonal validation datasets, and short-read-based transcriptomics to characterize placental transcript diversity and its relevance to metabolic and endocrine function in pregnancy.

---

## Placenta assembly annotation files

📄 [GTF of all transcripts assembled by ESPRESSO](Assembly/ESPRESSO_corrected_corrected.gtf.gz)

📄 [SQANTI3 classification table of all transcripts assembled by ESPRESSO](Assembly/ESPRESSO_corrected_classification.txt.gz)

📄 [GTF of filtered, high-confidence transcripts](Assembly/ESPRESSO_corrected_SC_filtered.gtf.gz)

📄 [SQANTI3 classification table of filtered, high-confidence transcripts](Assembly/ESPRESSO_corrected_SC_filtered_classification.txt.gz)

## Repository Structure

### I. Preparation of Datasets for Assembly  
📄 [`I_dataset_preparation.md`](I_dataset_preparation.md)  
Prepares CAGE-seq, DNase-seq, poly(A) site, and short-read splice junction data to validate long-read transcript structures.

### II. Transcriptome Assembly  
📄 [`II_assembly.md`](II_assembly.md)  
Processes Oxford Nanopore data through basecalling, trimming, correction, alignment, and de novo transcriptome assembly using ESPRESSO.

### III. Annotation, Quality Control, and Analysis  
📄 [`III_annotation_and_QC.md`](III_annotation_and_QC.md)

- 📜 [`IIIB_raw_assembly_analysis.R`](IIIB_raw_assembly_analysis.R)  
  Compares assembly to GENCODE v45 and GTEx v9 across gene/isoform counts and structural categories.

- 📜 [`IIIC_assembly_QC_and_filtering.R`](IIIC_assembly_QC_and_filtering.R)  
  Applies SQANTI3-based multi-layered filtering with CAGE-seq, poly(A) support, and splice junction evidence. Outputs filtered transcriptome and QC plots.

### IV. Short-Read Analyses in Birth Cohorts  
📄 [`IV_short_read_analyses.md`](IV_short_read_analyses.md)  
Explores the biological and clinical relevance of the filtered assembly using RNA-seq data from GUSTO (N=200) and Gen3G (N=152).

- 📜 [`IVB_short-read_txomics.R`](IVB_short-read_txomics.R)  
  Assesses quantification accuracy, expression patterns, and coding potential vs GENCODE.

- 📜 [`IVC_DGE_DTE.R`](IVC_DGE_DTE.R)  
  Identifies differentially expressed genes and transcripts associated with gestational diabetes mellitus (GDM).

- 📜 [`IVD_mediation_and_GOEA.R`](IVD_mediation_and_GOEA.R)  
  Tests for mediation of GDM effects on birth weight by transcript expression, and performs GO enrichment analysis.

---

## Key Outputs

- **Transcriptome Annotation**: Filtered GTF, FASTA, classification tables  
- **Figures**: All primary and supplementary figures (Figures 1–4, S1–S8)  
- **Expression Data**: Gene and transcript-level quantifications for GUSTO and Gen3G  
- **Mediation & Enrichment Results**: Individual and PC-based mediators, enriched GO terms

---

## Contact

For questions, please contact:  

📧 Sean T Bresnahan — [stbresnahan@mdanderson.org](stbresnahan@mdanderson.org)

🧬 [Bhattacharya Lab for Computational Genomics](https://bhattacharya-lab.com)

---
