#!/bin/bash
#SBATCH --job-name=trinotate_root
#SBATCH --output=trinotate_root_%j.out
#SBATCH --error=trinotate_root_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=standard


cd /anvil/scratch/x-kguse/root_analysis

module load biocontainers
module load trinotate

# Initialize the SQLite database (only run --create once)
#Trinotate Trinotate.sqlite.sqlite init \
 #--gene_trans_map trinity_out_dir.Trinity.fasta.gene_trans_map \
 #--transcript_fasta trinity_out_dir.Trinity.fasta \
 #--transdecoder_pep trinity_out_dir.Trinity.fasta.transdecoder.pep

# Load databases
Trinotate Trinotate.sqlite.sqlite LOAD_swissprot_blastx blastx.outfmt6
Trinotate Trinotate.sqlite.sqlite LOAD_swissprot_blastp blastp.outfmt6
Trinotate Trinotate.sqlite.sqlite LOAD_pfam pfam.domtblout

Trinotate Trinotate.sqlite.sqlite report > trinotate_annotation_report.tsv



