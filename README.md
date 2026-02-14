# Circular RNA-Seq Analysis Pipeline

Comprehensive circular RNA (circRNA) discovery, differential expression, structural characterization, and motif/miRNA enrichment analysis using:

* CIRCexplorer2
* STAR
* circRNAprofiler
* DESeq2

---

## 📌 Overview

This repository contains a complete end-to-end workflow for:

1. circRNA detection from RNA-Seq data
2. Differential expression analysis
3. Structural characterization of circRNAs
4. Motif enrichment analysis (RBP motifs)
5. miRNA binding site prediction

The workflow is demonstrated using dataset:

> **GSE235850** — RNA-Seq profiling of circular RNAs in colorectal cancer

---

# 🧬 Biological Background

Circular RNAs (circRNAs) are covalently closed RNA molecules formed via **back-splicing**, where a downstream 5′ splice donor joins an upstream 3′ splice acceptor.

Key properties:

* No 5′ or 3′ ends
* Resistant to exonucleases
* High stability
* Tissue-specific expression

### Functional roles

* miRNA sponges
* RBP scaffolding
* Protein coding (subset)
* Disease biomarkers (cancer, neurological, cardiovascular)

---
# ⚙️ Pipeline Overview

```
Raw FASTQ
   ↓
Quality Control (Trim Galore)
   ↓
Genome Alignment (STAR)
   ↓
circRNA Detection (CIRCexplorer2)
   ↓
Import into R (circRNAprofiler)
   ↓
Filtering
   ↓
Differential Expression (DESeq2)
   ↓
Structural Analysis
   ↓
Sequence Retrieval
   ↓
RBP Motif Enrichment
   ↓
miRNA Binding Prediction
```

---

# 🛠 Environment Setup

```bash
conda activate rnaseq_env
conda install -c bioconda circexplorer2 star trim-galore fastqc -y
```

R packages:

```r
BiocManager::install("circRNAprofiler")
BiocManager::install("DESeq2")
library(circRNAprofiler)
library(DESeq2)
```

---

# 🧹 Step 1: Quality Control

```bash
trim_galore \
  --fastqc \
  --paired \
  --cores 20 \
  sample_R1.fastq.gz \
  sample_R2.fastq.gz \
  -o Trimmed/
```

---

# 🧭 Step 2: Alignment with STAR

```bash
STAR \
  --genomeDir Genome_Index/ \
  --runThreadN 20 \
  --readFilesIn R1_val_1.fq.gz R2_val_2.fq.gz \
  --readFilesCommand zcat \
  --chimSegmentMin 10 \
  --chimJunctionOverhangMin 10 \
  --outSAMtype BAM SortedByCoordinate
```

**Important parameters:**

| Parameter                      | Purpose                            |
| ------------------------------ | ---------------------------------- |
| `--chimSegmentMin 10`          | Detect short back-splice junctions |
| `--chimJunctionOverhangMin 10` | Improves circRNA confidence        |
| `Chimeric.out.junction`        | Required for CIRCexplorer2         |

---

# 🔎 Step 3: circRNA Detection

```bash
fast_circ.py parse \
  -r hg38_kg.txt \
  -g hg38.fa \
  -t STAR \
  -o output_dir \
  Chimeric.out.junction
```

Output:

```
circularRNA_known.txt
```

---

# 📊 Downstream Analysis in R

## Initialize Project

```r
initCircRNAprofiler(
  projectFolderName = "projectCirc",
  detectionTools = "circexplorer2"
)
```

---

## Filter circRNAs

```r
filteredCirc <- filterCirc(mergedBSJunctions,
                           allSamples = FALSE,
                           min = 5)
```

---

# 📈 Differential Expression

Using DESeq2 via circRNAprofiler:

```r
deseqRes <- getDeseqRes(
    filteredCirc,
    condition = "A-B",
    fitType = "local",
    pAdjustMethod = "BH"
)
```

Volcano plot:

```r
volcanoPlot(deseqRes,
            log2FC = 1,
            padj = 0.05)
```

---

# 🧬 Structural Characterization

Functions used:

* `annotateBSJs()`
* `plotLenIntrons()`
* `plotLenBSEs()`
* `plotHostGenes()`
* `plotExBetweenBSEs()`
* `plotExPosition()`
* `plotTotExons()`

These compare:

* Flanking intron length
* Exon length
* Exon count
* Host gene complexity

---

# 🧬 Sequence Retrieval

Using BSgenome:

```r
library(BSgenome.Hsapiens.UCSC.hg19)
genome <- getBSgenome("BSgenome.Hsapiens.UCSC.hg19")
```

Retrieve:

| Function              | Purpose                |
| --------------------- | ---------------------- |
| `getCircSeqs()`       | Full circRNA sequences |
| `getSeqsAcrossBSJs()` | Junction sequences     |
| `getSeqsFromGRs()`    | Flanking regions       |

---

# 🧠 Motif Enrichment (RBP)

```r
motifs <- getMotifs(
    targets,
    width = 6,
    database = "ATtRACT",
    species = "Hsapiens",
    rbp = TRUE
)
```

Then:

```r
mergedMotifs <- mergeMotifs(motifs)
plotMotifs(...)
```

Normalization by sequence length is recommended.

---

# 🎯 miRNA Binding Site Prediction

```r
getMiRsites(targets)
```

Analyzes:

* Seed region (2–8)
* Central region (9–12)
* Compensatory region (13–16)

Supports:

* Canonical matches (A-U, G-C)
* G-U wobble
* Bulges/mismatches

---

# 📁 Expected Project Structure

```
projectCirc/
│
├── circexplorer2/
│   ├── circRNAs_001.txt
│   ├── circRNAs_002.txt
│   └── ...
│
├── experiments.txt
├── gencode.gtf
└── README.md
```

---

# 📊 Output

* Differentially expressed circRNAs
* Structural feature comparison plots
* RBP motif enrichment tables
* miRNA binding predictions
* Sequence FASTA files
