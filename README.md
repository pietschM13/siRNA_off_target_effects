# Computational analysis of miRNA-like off-target effects of RNA therapeutics

## Overview

This project investigates potential miRNA-like off-target effects of RNA therapeutics, with an initial focus on siRNA-mediated gene silencing.

siRNAs are designed to induce sequence-specific gene silencing through the RNA-induced silencing complex (RISC). However, the guide strand can also interact with partially complementary sequences in other transcripts, particularly through its seed region. These interactions may induce unintended repression of genes beyond the intended therapeutic target.

The project aims to investigate this phenomenon computationally by combining transcriptomic expression data with sequence-based analysis of potential seed-binding sites in annotated 3′UTRs.

The initial work focuses on one siRNA treatment and one RNA-seq dataset as a proof of concept. The analytical framework is intended to be scalable to additional siRNA and shRNA datasets.

## Research question

Do genes containing sequences complementary to the seed region of a therapeutic siRNA exhibit a systematically different expression distribution compared with genes without these sequences?

More specifically, is the presence of canonical miRNA-like seed-binding sites in annotated 3′UTRs associated with stronger gene repression following siRNA treatment?

## Objectives

1. **Dataset identification and acquisition**
   - Identify suitable RNA-seq datasets involving siRNA/shRNA treatments.
   - Retrieve and inspect expression data and associated experimental information.
   - Establish a reproducible data-processing workflow.

2. **On-target and differential expression analysis**
   - Assess the intended target gene's expression response.
   - Identify transcriptome-wide expression changes following treatment.
   - Estimate gene-level log2 fold changes using appropriate statistical methods.

3. **Sequence-based seed-site analysis**
   - Extract the guide-strand seed sequence from the therapeutic siRNA.
   - Retrieve annotated human 3′UTR sequences.
   - Identify canonical and non-canonical seed matches.
   - Assign seed-site presence/absence at gene level.

4. **Statistical analysis of potential off-target effects**
   - Compare log2 fold-change distributions between genes with and without seed matches.
   - Evaluate different seed-site definitions.
   - Assess whether seed-bearing genes show evidence of systematic repression.

5. **Visualisation and communication**
   - Visualise expression distributions and potential seed-mediated effects.
   - Document the analytical workflow and findings.
   - Explore scaling the analysis into an interactive web application.

## Initial study design

**Dataset:** GEO GSE164413

- Organism: *Homo sapiens*
- Cell model: Human adipose-derived stem cells
- Treatment: JARID2-targeting siRNA
- Control: Non-targeting siRNA control
- Experimental timepoints: Day 6 and day 13
- Replicates: Three matched control/treatment experiments per timepoint

The initial proof of concept focuses on day 6, representing an intermediate stage of differentiation.

**siRNA sequence:**

Reported passenger-strand sequence:

`GAAGAACGGGUGGUACGUA`

Candidate guide strand (reverse complement):

`UACGUACCACCCGUUCUUC`

The guide-strand assignment is provisional and based on sequence orientation rather than direct evidence of RISC loading.

**Initial seed definition:**

- Guide seed (positions 2–8): `ACGTACC`
- Canonical 8mer-A1 target-site motif: `AGGTACGT`

Additional site definitions, including 7mer-A1, 7mer-m8 and 6mer, will be considered subsequently.

## Data sources

**Gene expression**

- NCBI Gene Expression Omnibus (GEO)
- GSE164413 processed normalized and raw count matrices

**3′UTR sequences**

- TargetScan UTR sequence database
- Human records identified using NCBI Taxonomy ID 9606
- Database release/annotation approximately 2021

TargetScan sequences are retained at transcript level to account for alternative 3′UTR isoforms.

## Tools and technologies

- Python
- Pandas and NumPy
- Biopython
- SciPy and Statsmodels
- Matplotlib and Seaborn
- Jupyter Notebooks
- Visual Studio Code
- Linux / WSL
- Git and GitHub
- DESeq2 (planned for formal differential expression analysis)

## Analytical workflow

1. Dataset identification and acquisition.
2. Data inspection, cleaning and organisation.
3. On-target expression sanity check.
4. Guide-strand and seed-sequence definition.
5. Retrieval and processing of human 3′UTR sequences.
6. Transcript-level seed-site identification.
7. Aggregation of seed-site annotations to gene level.
8. Integration with differential expression results.
9. Statistical comparison of log2 fold-change distributions.
10. Visualisation, interpretation and documentation.

## Statistical approach

The central analysis will compare the cumulative distributions of gene-level log2 fold changes between:

- Genes containing a particular seed-binding site.
- Genes without that seed-binding site.

Wilcoxon rank-sum testing will be used to assess differences between groups, with complementary distribution-based analyses considered where appropriate.

The analysis will include genes with measurable expression and usable 3′UTR annotations, rather than restricting the population to significantly differentially expressed genes.

Potential confounding factors, including 3′UTR length, GC content, expression level and transcript isoform diversity, will be considered during interpretation and subsequent model development.

Importantly, a sequence match represents a potential targeting interaction, not proof of direct binding or functional repression.

## Current scope and future development

The initial objective is to establish a complete analytical proof of concept using one siRNA treatment and one dataset.

Subsequent development may include:

- Additional siRNA and shRNA datasets.
- Multiple therapeutic sequences and seed definitions.
- Cross-dataset comparisons.
- More advanced statistical modelling.
- An interactive application for displaying and exploring potential off-target effects.
