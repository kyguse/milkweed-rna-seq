# RNA-Seq Pipeline Script

This folder contains a general-purpose pipeline script used for the RNA-Seq analysis of milkweed samples (leaf, root, and hybrid).

## 📄 File

- `rna_seq_pipeline_milkweed.sh`: A step-by-step Bash-style script outlining the key stages of the RNA-Seq workflow from raw reads to transcript quantification and downstream analysis.

## 🧠 Intended Use

This script is designed for **teaching and reproducibility**. It highlights major steps used in de novo transcriptome analysis and annotation, including:

- Quality control with FastQC
- Transcriptome assembly using Trinity
- Read alignment with Bowtie2
- Quantification with Salmon
- Functional annotation with BLAST and Trinotate
- Differential expression analysis with DESeq2 and GO enrichment using TopGO

Each step is documented with commands and brief comments. Students can adapt these commands for their own data and environments.

## 📚 Tutorial Resource

For a guided video walkthrough of this pipeline, see:

▶️ [RNA-seq Pipeline Overview (YouTube)](https://youtu.be/q8iTNdVWpLQ?si=CtPJ8W66ljIZekI7)

## 📁 How to Use

Clone this repository or download the script directly, and execute on a high-performance computing environment like ACCESS Anvil. Modules and file paths should be adjusted to match your local setup.

```bash
bash rna_seq_pipeline_milkweed.sh
```

---

Created by Kylene Guse, 2025
