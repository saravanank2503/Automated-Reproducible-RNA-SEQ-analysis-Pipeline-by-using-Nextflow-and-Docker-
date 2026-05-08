# Automated-Reproducible RNA-Seq Analysis Pipeline

![Pipeline](https://img.shields.io/badge/Pipeline-RNA--Seq-blue) ![Version](https://img.shields.io/badge/Version-2.0.0-green) ![Platform](https://img.shields.io/badge/Platform-macOS%20%7C%20Linux-lightgrey) ![Docker](https://img.shields.io/badge/Docker-Required-blue)

## What Is This?

A fully automated, containerized RNA-Seq analysis pipeline built with **Nextflow** and **Docker**. It takes raw paired-end FASTQ files and produces publication-ready differential expression results — from quality control through to DESeq2 statistical analysis and visualization — with a single command.

## Why It Matters

Understanding which genes are turned on or off across conditions is central to modern biology. This pipeline enables:

- **Reproducible science** — every tool is version-pinned in a Docker container; results are identical across machines
- **End-to-end automation** — no manual steps between raw reads and DEG results
- **Clinical-grade reference** — uses GRCh38 / GENCODE v44, the current human genome standard
- **Scalable design** — handles any number of samples and conditions with automatic pairwise comparisons

This pipeline was validated on a **CAR-T cell immunotherapy dataset** (4 conditions × 3 replicates), comparing OCI-AML3, THP-1, CAR-T, and control samples to identify transcriptomic signatures relevant to anti-tumor response.

---

## Pipeline Overview

```
FASTQ → FastQC → Trim Galore → STAR Align → Samtools → featureCounts → Merge → DESeq2 → MultiQC
```

| Step | Tool | Version |
|------|------|---------|
| Raw QC | FastQC | 0.12.1 |
| Trimming | Trim Galore | 0.6.10 |
| Alignment | STAR | 2.7.11b |
| BAM sort/index | Samtools | 1.19 |
| Counting | featureCounts | 2.0.1 |
| Differential expression | DESeq2 | Bioconductor 3.19 |
| QC Report | MultiQC | 1.21 |

---

## System Requirements

| Resource | Minimum | Recommended |
|----------|---------|-------------|
| RAM | 32 GB | 48 GB |
| CPU | 4 cores | 6–8 cores |
| Disk | 200 GB | 500 GB |
| Docker | 24.0+ | latest |
| Nextflow | 23.04+ | latest |
| Java | 11+ | 17+ |

---

## Project Structure

```
rnaseq_pipeline/
├── main.nf                    ← Nextflow workflow
├── nextflow.config            ← Resources, containers, profiles
├── samplesheet.csv            ← Sample metadata
├── data/                      ← Raw FASTQ files
├── modules/                   ← One process per file
│   ├── fastqc.nf
│   ├── trimgalore.nf
│   ├── star_align.nf
│   ├── samtools.nf
│   ├── featurecounts.nf
│   ├── merge_counts.nf
│   ├── deseq2.nf
│   └── multiqc.nf
├── bin/                       ← Standalone R scripts
│   ├── deseq2_analysis.R
│   └── merge_counts.R
├── docker/                    ← Custom Docker build
│   └── Dockerfile
└── results/                   ← Pipeline outputs
```

---

## Quick Start

### 1. Set Environment Variables

```bash
export PROJECT_DIR=/path/to/your/rnaseq_project
export REF_DIR=/path/to/your/reference
```

### 2. Pull Docker Images

```bash
docker pull quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0
docker pull quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0
docker pull quay.io/biocontainers/star:2.7.11b--h43eeafb_0
docker pull quay.io/biocontainers/samtools:1.19--h50ea8bc_0
docker pull quay.io/biocontainers/subread:2.0.1--hed695b0_0
docker pull quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0
```

### 3. Build the DESeq2 Container

```bash
cd $PROJECT_DIR/docker
docker build --platform linux/amd64 --no-cache -t rnaseq_deseq2:1.0 .
```

### 4. Prepare Reference Files

Reference files should be placed at `$REF_DIR/`:
```
$REF_DIR/
├── genome.fa                       ← GRCh38 primary assembly FASTA
├── star_index/                     ← Pre-built STAR index (sjdbOverhang 149)
└── gtf/
    └── gencode.v44.annotation.gtf  ← Gene annotation
```

Download if needed:
```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
gunzip *.gz
```

> **Important:** FASTA and GTF must be from the same source and version. Mixing sources (e.g. GENCODE FASTA + Ensembl GTF) causes featureCounts assignment failures.

### 5. Prepare Your Samplesheet

```csv
sample,fastq_1,fastq_2,condition,replicate,batch
CTRL_1,data/SAMPLE_1_R1.fastq.gz,data/SAMPLE_1_R2.fastq.gz,control,1,batch1
TREAT_1,data/SAMPLE_2_R1.fastq.gz,data/SAMPLE_2_R2.fastq.gz,treatment,1,batch1
```

- Minimum **3 replicates per condition** for DESeq2
- Sample names must be unique, no spaces
- Condition names are case-sensitive

### 6. Run the Pipeline

```bash
cd $PROJECT_DIR

nextflow run main.nf \
    -profile docker \
    --samplesheet samplesheet.csv \
    --star_index  $REF_DIR/star_index \
    --gtf         $REF_DIR/gtf/gencode.v44.annotation.gtf \
    --outdir      results \
    --max_cpus    6 \
    --max_memory  48.GB
```
---

## Expected Outputs

```
results/
├── fastqc/                         ← Raw and trimmed QC reports
├── trimgalore/                     ← Trimmed FASTQs and logs
├── bam/                            ← Sorted, indexed BAM files
├── counts/
│   └── merged_counts.txt           ← Count matrix (input to DESeq2)
├── deseq2/
│   ├── deseq2_results_*.csv        ← DEG results per comparison
│   ├── deseq2_normalized_counts.csv
│   └── plots/
│       ├── pca.pdf
│       ├── dispersion.pdf
│       ├── sample_distance_heatmap.pdf
│       ├── top50_degs_heatmap.pdf
│       ├── volcano_*.pdf           ← One per pairwise comparison
│       └── maplot_*.pdf
└── multiqc/
    └── multiqc_report.html         ← Start here for QC overview
```

### Quality Thresholds

| Step | Pass | Warning | Fail |
|------|------|---------|------|
| Trimming | >80% reads kept | 60–80% | <60% |
| STAR mapping | >70% unique | 50–70% | <50% |
| featureCounts assigned | >50% | 30–50% | <30% |

---

## Key Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `min_length` | 36 bp | Min read length after trimming |
| `quality_cutoff` | 20 | Phred quality threshold |
| `adapter` | `AGATCGGAAGAGC` | Illumina TruSeq — always set explicitly |
| `strandedness` | 2 | 0 = unstranded, 1 = stranded, 2 = reverse |
| `deseq2_fdr` | 0.05 | FDR threshold |
| `deseq2_lfc` | 1.0 | log2 fold-change threshold |
| `max_memory` | 48 GB | Total Docker RAM allocation |
| `max_cpus` | 6 | CPU cores |

---

## Troubleshooting

### Find a failed task
```bash
nextflow log <run_name> -f name,status,exit | grep FAILED
cat work/xx/xxxxxxxx/.command.err
```

### Common exit codes

| Code | Cause | Fix |
|------|-------|-----|
| 125 | Docker crashed or OOM | Restart Docker; verify RAM ≥ 40 GB |
| 137 | Process killed by OS | Increase Docker memory limit |
| 127 | Command not found in script | Check for broken line continuations in `modules/*.nf` |
| 2 | Low mapping / wrong adapter | Verify adapter sequence and strandedness setting |
| 1 | R/DESeq2 error | Check container; verify sample names match samplesheet |

### When NOT to use `-resume`
If you changed the adapter, genome, GTF, or replaced input files — delete `.nextflow/cache` and `work/` and do a fresh run.

---

## Design Principles

```
Layer 1: Nextflow   → orchestration (what runs, in what order)
Layer 2: Docker     → execution environment (reproducibility)
Layer 3: R scripts  → statistical analysis (bin/*.R)
```

- All containers are version-pinned — never use `:latest`
- R logic lives in `bin/` only, never inside Nextflow script blocks
- STAR runs with `maxForks = 1` —  this ensure at the time 1 sample only process for alignment, if your RAM is high you can do simultaneously alignment for samples
---
*Automated-Reproducible Pipeline  — Nextflow + Docker + DESeq2 | GRCh38 / GENCODE v44*
