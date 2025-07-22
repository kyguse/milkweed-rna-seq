#!/bin/bash
#SBATCH --job-name=hmmscan
#SBATCH --output=hmmscan_%j.out
#SBATCH --error=hmmscan_%j.err
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=standard

cd /anvil/scratch/x-kguse/root_analysis

module load hmmer/3.3.2


hmmscan --cpu 16  --domtblout pfam.domtblout /anvil/scratch/x-kguse/root_analysis/Pfam-A.hmm /anvil/scratch/x-kguse/root_analysis/trinity_out_dir.Trinity.fasta.transdecoder_dir/longest_orfs.pep
