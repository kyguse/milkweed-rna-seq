#!/bin/bash
#SBATCH -t 06:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --job-name=trinotate
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --error=%x-%J-%u.err
#SBATCH --output=%x-%J-%u.out


ml biocontainers trinotate/3.2.2


Trinotate /anvil/scratch/x-kguse/leaf/Trinotate.sqlite.sqlite init \
  --gene_trans_map /anvil/scratch/x-kguse/leaf/trinity_out_dir.Trinity.fasta.gene_trans_map \
  --transcript_fasta /anvil/scratch/x-kguse/leaf/trinity_out_dir.Trinity.fasta \
  --transdecoder_pep /anvil/scratch/x-kguse/leaf/trinity_out_dir.Trinity.fasta.transdecoder.pep

Trinotate Trinotate.sqlite.sqlite LOAD_swissprot_blastx /anvil/scratch/x-kguse/leaf/blastx.outfmt6
Trinotate Trinotate.sqlite.sqlite LOAD_swissprot_blastp /anvil/scratch/x-kguse/leaf/blastp.outfmt6
Trinotate Trinotate.sqlite.sqlite LOAD_pfam /anvil/scratch/x-kguse/leaf/output.domtblout
Trinotate Trinotate.sqlite.sqlite report > trinotate_annotation_report.tsv




