# This will create a CSV file with the paths to the TB fasta files

#!/bin/bash

# Directory where GPSCs are stored
GPSC_dir="/vast/palmer/scratch/turner/flg9/snakemake_workflows/pangenome_alignment/GitHub/pipelines/prokka_roary/seqs/TB/"

# Output CSV file
output_file="TB_seqs.tsv"

# Write the header line
echo -e "sample_id\tseq_path" > $output_file

# Loop over the fasta files in GPSCs dir
for fasta_file in $GPSC_dir*.fa; do
    # Extract the sample_id from the file name
    sample_id=$(basename $fasta_file .fa)
    # Write a line to the CSV file with the path to the fasta file
    echo -e "$sample_id\t$fasta_file" >> $output_file
done