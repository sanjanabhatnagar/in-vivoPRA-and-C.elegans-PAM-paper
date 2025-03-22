# Name of the Paper
This repository contains the scripts and the steps to execute these to extract 

# Command for running mafft on fasta files containing ortholog groups.
# Performing Percent identity check to keep 1:1 ortholog.

nohup bash run_multi_mafft.sh CDS_SwitchEvents ./Neuro_Mus_Switch/ > mafft_CDS.out &

#It runs mafft on all fasta files in parallel

~/miniconda3/envs/sbenv/bin/mafft --adjustdirection --quiet --auto --thread 3  $file > "${output_folder}/${directory}_aln/${filename}.aln"

#Two options generate reverse complement sequences, as necessary, 
#and align them together with the remaining sequences.
# is based on 6 mer counting and faster.

# The former works well in most cases, unless the sequences are highly diverged.
