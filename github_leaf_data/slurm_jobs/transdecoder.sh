#!/bin/bash
#SBATCH -t 06:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --job-name=transdecoder
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --error=%x-%J-%u.err
#SBATCH --output=%x-%J-%u.out

ml biocontainers transdecoder blast/2.13.0 hmmer/3.3.2
cd /anvil/scratch/x-kguse/leaf

hmmscan --domtblout output.domtblout Pfam-A.hmm longest_orfs.pep

TransDecoder.Predict -t /anvil/scratch/x-kguse/leaf trinity_out_dir.Trinity.fasta \
	--retain_pfam_hits output.domtblout \
	--retain_blastp_hits blastp.outfmt6



