#!/bin/bash

# point to dir where seqs are stored
TB_GPSC_dir="/vast/palmer/scratch/turner/flg9/snakemake_workflows/pangenome_alignment/GitHub/pipelines/genome_annotation/seqs/TB/"
SP_GPSC_dir="/vast/palmer/scratch/turner/flg9/snakemake_workflows/pangenome_alignment/GitHub/pipelines/genome_annotation/seqs/SP/"

# specify output file
output_file="SP_TB_seqs.tsv"

# create header for output file
echo -e "sample_id\tseq_path\tgenus" > $output_file

# loop over TB fasta files and add "Mycobacterium" as genus
for fasta_file in $TB_GPSC_dir*.fa; do
    sample_id=$(basename $fasta_file .fa)
    echo -e "$sample_id\t$fasta_file\tMycobacterium" >> $output_file
done

# Loop over SP fasta files and add "Streptococcus" as genus
for fasta_file in $SP_GPSC_dir*.fasta; do
    sample_id=$(basename $fasta_file .fasta)
    echo -e "$sample_id\t$fasta_file\tStreptococcus" >> $output_file
done