#!/bin/bash
#SBATCH --job-name=blastp
#SBATCH --output=blastp_%j.out
#SBATCH --error=blastp_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G


#change directory
cd  /anvil/scratch/x-kguse/root_analysis

# Load BLAST module
module load biocontainers
module load blast/2.13.0

# run blastp
blastp -query /anvil/scratch/x-kguse/root_analysis/trinity_out_dir.Trinity.fasta.transdecoder_dir/longest_orfs.pep \
    -db uniprot_sprot.fasta  -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads 10 > blastp.outfmt6
