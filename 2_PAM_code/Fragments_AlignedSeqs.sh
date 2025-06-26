#!/bin/bash

tr_name_file=$1  #This file contains transcript names, sequences of which are used for the alignments.Since, we only want to use coordinates as per these transcripts from the coordinates file generated using .gff.
directory=$2  # Directory with longest transcript alignments.
output_folder=$3

echo "Directory for search: $directory"
echo "dir ${output_folder}/${directory}"
mkdir -p ${output_folder}/introns_${directory}

while IFS= read -r tr_name; do
    echo "Searching aligned files for : $tr_name"

    for file in "$output_folder/$directory"/*.aln; do
            filename=$(basename "$file" | awk -F "_PIDfilter.fa.aln" '{print $1}')
        if [[ "$tr_name" == *"$filename"* ]]; then
                echo "The $tr_name matched with $filename fasta :o. Now extracting intron fragments! :)"
             python IntronSegments_alnfastas.py "$tr_name" "$file" ./Exon_IntronFragment_coordinates.tsv "${output_folder}/introns_${directory}/"
        fi
    done

done < "$tr_name_file"
