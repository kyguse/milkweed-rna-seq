# 🧬 Milkweed RNA-Seq Pipeline

**De novo transcriptome assembly and differential expression analysis of *Asclepias* species and their hybrid**

---

## 🌿 Project Overview

This project analyzes root and leaf transcriptomes from *Asclepias speciosa*, *A. syriaca*, and their F1 hybrid as part of a study on microevolutionary responses to environmental variation. The analysis includes quality control, de novo assembly, transcript quantification, differential expression testing, and gene ontology enrichment.

> 💡 Conducted on ACCESS Anvil HPC  
> 🔬 Tools: Trinity, Bowtie2, Salmon, Transdecoder, Trinity, DESeq2, TopGO  
> 📁 Data: Paired-end RNA-seq (Illumina)

## 🚧 Project Status: In Progress
This analysis is currently ongoing. Documentation, results, and scripts are being actively developed and refined.


---

## 🛠 Pipeline Overview

```bash
FASTQ → Trimmomatic → Trinity → Bowtie2 → ncbi-blast → Salmon → Transdecoder → Trinotate → in R →  DESeq2 → TopGO


---

## 🛠 Tools & Versions

| Tool        | Version   | Description                           |
|-------------|-----------|---------------------------------------|
| Trinity      | v2.15.1   | De novo transcriptome assembly        |
| Salmon       | v1.10.1   | Transcript quantification             |
| Bowtie2      | v2.5.1    | Read alignment (optional)             |
| DESeq2       | R package | Differential expression analysis      |
| TopGO        | R package | Gene Ontology enrichment              |
| Transdecoder | v5.5.0    | Identify candidate coding regions     |
| Trinotate    | v3.2.2    | Functional annotation                 |

---

## 📂 Directory Structure

milkweed-rna-seq/
├── data/ # sample metadata and future path raw fastq files
├── scripts/ # Bash and R scripts
├── slurm_jobs/ # SLURM job submission files
├── trinity_out/ # Trinity assemblies
├── quants/ # Salmon quant outputs
├── results/ # DESeq2 outputs, plots, tables
└── README.md # This file



