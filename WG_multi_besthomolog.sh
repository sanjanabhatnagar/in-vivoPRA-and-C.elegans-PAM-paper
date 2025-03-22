#!/bin/bash

directory=$1  # Directory with Best homologs aa alignments
directory2=$2 # Directory with Whole gene sequences
output_folder=$3

echo "Directory for search: $directory2"
echo "Creating output directory: ${output_folder}/BH_${directory2}"
mkdir -p "${output_folder}/BH_${directory2}"

for file in "$output_folder/$directory"/*best_homologs.fa; do
    filename=$(basename "$file" | awk -F "_best_homologs.fa" '{print $1}')
    echo "------------Working on file $filename --------------"
    for file2 in "$output_folder/$directory2"/*.fa; do
            wg_filename=$(basename "$file2" | awk -F ".fa" '{print $1}')
            echo "++++++++++++++$wg_filename+++++++++++++++"
            if [[ "$filename" == *"$wg_filename"* ]]; then
                    echo "Finding best homologs in $wg_filename fasta based on $filename amino acid alignments PID."
                    python WG_besthomolog.py "$file" "$file2" "${output_folder}/BH_${directory2}/${wg_filename}_PIDfilter.fa"
            fi
    done
done
