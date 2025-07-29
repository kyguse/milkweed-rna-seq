#!/bin/bash
#SBATCH -t 06:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --job-name=trinotate
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --error=%x-%J-%u.err
#SBATCH --output=%x-%J-%u.out

ml biocontainers
ml blast

blastx -db unisprot.pep\
	query trinity_out_dir.Trinity.fasta -num_threads 2 \
	max_target_seqs 1 -outfmt 6 -evalue 1e-5 \
	> swissprot.blastx.outfmt6

