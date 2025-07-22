#!/bin/bash

# ==========================================
# RNA-Seq Analysis Pipeline for Milkweed Project
# Author: Kylene Guse
# Description: This script summarizes the RNA-seq pipeline used for leaf, root, and hybrid samples
# ==========================================

# STEP 1: Run FastQC on raw FASTQ files
# Example:
# fastqc sample_R1.fastq.gz sample_R2.fastq.gz

# STEP 2: Concatenate left and right reads for Trinity
# cat trimmed_R1/* > all_R1.fastq.gz
# cat trimmed_R2/* > all_R2.fastq.gz

# STEP 3: Load Trinity and Samtools modules
# module load trinity/2.15.1
# module load samtools

# STEP 4: Run Trinity on the concatenated paired-end reads
# Trinity --seqType fq --left all_R1.fastq.gz --right all_R2.fastq.gz --CPU 16 --max_memory 100G

# STEP 5: Assess Assembly Quality with TrinityStats
# TrinityStats.pl trinity_out_dir/Trinity.fasta

# STEP 6: Align reads to the assembly using Bowtie2
# bowtie2-build Trinity.fasta Trinity_index
# bowtie2 -x Trinity_index -1 all_R1.fastq.gz -2 all_R2.fastq.gz | samtools view -Sb - > aligned.bam

# STEP 7: Functional Annotation with BLAST and Trinotate
# wget uniprot_sprot.fasta.gz
# gunzip uniprot_sprot.fasta.gz
# makeblastdb -in uniprot_sprot.fasta -dbtype prot
# blastx -query Trinity.fasta -db uniprot_sprot.fasta -out blastx.outfmt6 -evalue 1e-20 -outfmt 6 -num_threads 8

# STEP 8: Quantify transcripts using Salmon
# salmon index -t Trinity.fasta -i salmon_index
# salmon quant -i salmon_index -l A -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -p 8 -o quants/sample

# STEP 9: Prepare samples file for tximport
# Format: condition<TAB>sample_id

# STEP 10: Differential Expression and GO Analysis
# Run DESeq2 and TopGO in R

# Additional Notes:
# See video tutorial for pipeline explanation:
# https://youtu.be/q8iTNdVWpLQ?si=CtPJ8W66ljIZekI7
