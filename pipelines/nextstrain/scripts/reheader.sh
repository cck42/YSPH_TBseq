#!/bin/bash

# Input file containing the list of names
input_file="data/seqnames.txt"

# Output TSV file
output_file="output.tsv"

# Create or clear the output file
> "$output_file"

# Read the input file line by line
while IFS= read -r line; do
    # Extract the shortened name (everything after 'results/' and before the second '/')
    shortened_name=$(echo "$line" | awk -F'/' '{print $3}')
    shortened_name_2=$(echo "$shortened_name" | awk -F'_' '{print $1}')
    
    # Write the old name and the shortened name to the output file in TSV format
    echo -e "$line\t$shortened_name_2" >> "$output_file"
done < "$input_file"

echo "Output written to $output_file"