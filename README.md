# 🧬 Milkweed RNA-Seq Pipeline

**De novo transcriptome assembly and differential expression analysis of *Asclepias* species and their hybrid**

---

## 🌿 Project Overview

This project analyzes root and leaf transcriptomes from *Asclepias speciosa*, *A. syriaca*, and their F1 hybrid as part of a study on microevolutionary responses to environmental variation. The analysis includes quality control, de novo assembly, transcript quantification, differential expression testing, and gene ontology enrichment.

> 💡 Conducted on ACCESS Anvil HPC  
> 🔬 Tools: Trinity, Bowtie2, Salmon, Transdecoder, Trinity, DESeq2, TopGO  
> 📁 Data: Paired-end RNA-seq (Illumina)

---

## 🛠 Pipeline Overview

```bash
FASTQ → Trimmomatic → Trinity → Bowtie2 → ncbi-blast → Salmon → Transdecoder → Trinotate → in R →  DESeq2 → TopGO
