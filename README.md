# Name of the Paper


## Performing protein alignments to find best orthologs

The protein alignments were conducted using MAFFT (Katoh et al., 2019). MAFFT can be installed using conda - 
```conda install -c bioconda mafft
```
MAFFT was used with the following parameters. This can be enclosed within shell scripts to parallely run on multiple protein sequence files. 

```
~/mafft --quiet --auto --thread 3  protein_seq > protein_seq.aln
```
Next, Best_homologfinder.py is ran on these protein sequence alignments, and this picks the best ortholog per species for C. elegans from the ortholog group. This returns a final output file with the IDs of best orthologs that can be further used to filter the nucleotide seqeuence fasta files.



#It runs mafft on all fasta files in parallel

~/miniconda3/envs/sbenv/bin/mafft --adjustdirection --quiet --auto --thread 3  $file > "${output_folder}/${directory}_aln/${filename}.aln"

#Two options generate reverse complement sequences, as necessary, 
#and align them together with the remaining sequences.
# is based on 6 mer counting and faster.

# The former works well in most cases, unless the sequences are highly diverged.

# References

Katoh, K., Rozewicki, J., & Yamada, K. D. (2019). MAFFT online service: multiple sequence alignment, interactive sequence choice and visualization. 20(4), 1160–1166. https://doi.org/10.1093/bib/bbx108

Alam, A., Duncan, A. G., Mitchell, J. A., & Moses, A. M. (biorxiv). Functional similarity of non-coding regions is revealed in phylogenetic average motif score representations. https://doi.org/10.1101/2023.04.09.536185


