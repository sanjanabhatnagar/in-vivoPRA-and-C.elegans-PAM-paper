#!/bin/bash

directory=$1  # Directory with amino acid alignments
output_folder=$2

echo "Directory for search: $directory"
echo "dir ${output_folder}/${directory}"
mkdir -p "${output_folder}/BH_${directory}"

for file in "$output_folder/$directory"/*.aln; do
        filename=$(basename "$file" | awk -F ".fa.aln" '{print $1}')
        echo "Finding best homologs in  $filename fasta."
        python Best_homologfinder.py "$file"  "${output_folder}/BH_${directory}/${filename}_best_homologs.fa"
done
