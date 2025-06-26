#!/bin/bash

genes=$1  #This file contains transcript names on either +strand or on -ve strand.
strand=$2
directory=$3  # Directory with extracted Intron fragments.
output_folder=$4

echo "Directory for reverse_complement: $directory"
echo "dir ${output_folder}/${directory}"

while IFS= read -r Strand; do
    echo "Searching aligned files for : $Strand"
    for file in "$output_folder/$directory"/*.fasta; do
        filename=$(basename "$file" | cut -f 1 -d '.')
        if [[ "$Strand" == *"$filename"* ]]; then
             python revcomp_afterIntronFragments.py "$file" "$strand" "$file.reverse_comp"
        fi
    done

done < "$genes"              
