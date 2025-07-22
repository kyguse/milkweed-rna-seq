#!/bin/bash
#SBATCH --job-name=transdecoder
#SBATCH --output=transdecoder_%j.out
#SBATCH --error=transdecoder_%j.err
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=standard

cd /anvil/scratch/x-kguse/root_analysis
module load biocontainers
module load transdecoder

TransDecoder.Predict -t /anvil/scratch/x-kguse/root_analysis/trinity_out_dir.Trinity.fasta \
--retain_pfam_hits pfam.domtblout \
--retain_blastp_hits blastp.outfmt6




