# Name of the Paper


## 1. Performing protein alignments to find best orthologs

The protein alignments were conducted using MAFFT (Katoh et al., 2019). MAFFT can be installed using conda - 
```
conda install -c bioconda mafft
```
MAFFT was used with the following parameters. This can be enclosed within shell scripts to parallely run on multiple protein sequence files. 

```
~/mafft --quiet --auto --thread 3  protein_seq > protein_seq.aln
```
Next, best_multi_homologfinder.sh is ran on these protein sequence alignments, and this picks the best ortholog per species for C. elegans from the ortholog group by calculating percent identity (PID). It calls the Best_homologfinder.py script and outputs - 

1. Best Ortholog, with highest PID  sequence containing files into the user specified directory ("*_best_homologs.fa" extension). For example - 
```
>CELEG.F36H1.2e
---------------------------------MV---------QT--------------
--LRR---------------P--W----------QEAA--SSAFAVASALPV-TMNSTQI
AELFEQVEHGTTE-LRCALTAEISALRNANGESLLTVAVRSGNTAVAKQLAQLDPDA-ID
...
>CAFRA.g8544.t2_PID:0.1400137268359643
---------------------------------MVV--------QT--------------
--LLR---------------P--W----------QEAA--TSAFAVASALPV-TMNSTQI
AELFEQVEQGESEQLRCALTAELISMRNANGESLLVVAARVGNSAVAKQLIHLESSQFLN
...
```
2. File containing list of all orthologs from a species and corresponding PID with C.elegans ortholog within the group.
```
Event	Ref	Sps	PID
./AminoAcid_SwitchEvents_aln/B0348.4.fa.aln	CELEG	CBECE.CSP29.g10568.t2	0.10133333333333333
./AminoAcid_SwitchEvents_aln/B0348.4.fa.aln	CELEG	CBOVI.g1040.t1	0.20364238410596028
./AminoAcid_SwitchEvents_aln/B0348.4.fa.aln	CELEG	CBOVI.g1041.t2	0.19103773584905662
./AminoAcid_SwitchEvents_aln/B0348.4.fa.aln	CELEG	CBREN.CBN18509	0.07034372501998401
./AminoAcid_SwitchEvents_aln/B0348.4.fa.aln	CELEG	CBREN.CBN25584	0.07344632768361582
...
```
## 2. Extracting best orthologs nucleotide sequences based on previously calculated protein alignment PIDs.

Next, WG_multi_besthomolog.sh is run to subset the whole gene nucleotide sequences and only include ortholog sequences with highest PID for a given species, per protein alignment. The shell script runs  WG_besthomolog.py on each file in the given directory.
```
nohup bash WG_multi_besthomolog.sh BH_protien_alignments_directory nucleotide_fasta_directory ./output_folder/ > output.txt &
```
It outputs files with *_PIDfilter.fa!

Additional steps are performed to ensure there's only one occurence of each species name in each file in the directory and thus only one, best ortholog is chosen from each species. 

## 3. Whole gene alignments within ortholog groups

The obtained gene sequence files with best orthologs, are aligned using MAFFT. The purpose behind aligning the whole gene sequences was to further extract introns in each species based on exon alignments. This time --adjustdirection parameter is used so that MAFFT can find the best orientation for a given orthologous sequence and align it to C.elegans ortholog per ortholog group.  Following command is used - 
```
~/mafft --adjustdirection --quiet --auto --thread 3  Gene_PIDfilter_file > Gene_PIDfilter_file.aln
```
Before proceeding the orientation of C.elegans sequence is checked and it is ensured that negative strand genes and positive strand genes are in correct orientation. 

## 4. Extracting per species intron fragments from whole gene sequence alignments.





# References

Katoh, K., Rozewicki, J., & Yamada, K. D. (2019). MAFFT online service: multiple sequence alignment, interactive sequence choice and visualization. 20(4), 1160–1166. https://doi.org/10.1093/bib/bbx108

Alam, A., Duncan, A. G., Mitchell, J. A., & Moses, A. M. (biorxiv). Functional similarity of non-coding regions is revealed in phylogenetic average motif score representations. https://doi.org/10.1101/2023.04.09.536185


