# Name of the Paper

# Table of contents - 


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

Next, WG_multi_besthomolog.sh is run to subset the longest transcript nucleotide sequences and only include ortholog sequences with highest PID for a given species, per protein alignment. The shell script runs  WG_besthomolog.py on each file in the given directory.
```
nohup bash WG_multi_besthomolog.sh BH_protien_alignments_directory nucleotide_fasta_directory ./output_folder/ > output.txt &
```
It outputs files with *_PIDfilter.fa!

Additional steps are performed to ensure there's only one occurence of each species name in each file in the directory and thus only the best ortholog is chosen from each species. 

## 3. Longest transcript alignments within ortholog groups

The obtained gene sequence files with best orthologs, are aligned using MAFFT. The purpose behind aligning the longest transcript sequences was to further extract introns in each species based on exon alignments. This time --adjustdirection parameter is used so that MAFFT can find the best orientation for a given orthologous sequence and align it to C.elegans ortholog per ortholog group. Additionally, we ensured that in each ortholog group fasta file, C.elegans sequence was the first as MAFFT forces the orientation of the first sequence of the sequence file it is aligning. Following command is used - 
```
~/mafft --adjustdirection --quiet --auto --thread 3  Gene_PIDfilter_file > Gene_PIDfilter_file.aln
```
Before proceeding the orientation of C.elegans sequence is checked and it is ensured that negative strand genes and positive strand genes are in correct orientation. Since, the exons and introns are located using coordinates, negative strand gene transcripts should be on negative strand whereas positive strand genes shoud be on positive strand (following the reference species, C.elegans) for the next step. 

## 4. Extracting per species intron fragments from longest transcript sequence alignments.

```
nohup bash Fragments_AlignedSeqs.sh ./Neuro_Mus_Switch/TranscriptNames_Cutterlab.txt BH_nucleotide_fasta_aln_directory ./output_folder/ > output.txt &

#This in turn calls IntronSegments_alnfastas.py and it is run on different ortholog groups as follows - 
#python IntronSegments_alnfastas.py B0348.4d.1 ./Neuro_Mus_Switch//BH_nucleotide_fasta_aln_directory/B0348.4_PIDfilter.fa.aln #~/Exon_IntronFragment_coordinates.tsv ~/introns_BH_nucleotide_fasta_aln_directory/
```
It outputs this file while running and the information can be checked to ensure it is running as expected.

```
Directory for search: BH_nucleotide_fasta_aln_directory
dir ./BH_nucleotide_fasta_aln_directory/
Searching aligned files for : B0348.4d.1
The B0348.4d.1 matched with B0348.4 fasta :o. Now extracting intron fragments! :)
The exon start : 0. The exon end : 83
       start  end
10424     84  154
The length of C.elegans intron fragment is : 70
The relative ungapped coordinates are 84 and 154
The corresponding gapped coordinates are 386 and 1751
*-------------------------------------*-*-------------------------------------*
The exon start : 1488. The exon end : 1849
       start   end
10425   1417  1487
10427   1850  1920
The length of C.elegans intron fragment is : 70
The relative ungapped coordinates are 1417 and 1487
The corresponding gapped coordinates are 6224 and 6464
...
```
The inference of adjacent introns relies on a given exon. Therefore, if an exon in one species cannot be aligned to exons in others, the adjacent introns for that exon will not be inferred in those species. The inferred introns are checked for fragments of same size across each species in the files chosen arbitrarily. In this case, we are using all introns shorter than 70 bp and the longer introns are broken into 70 bp fragments from either splice site. 

For calculating the conservation scores, we ensure that sequences within an ortholog group does not contain ambiguous nucleotides, doesn't contain redundant sequences and contains sequences with a minimum distance of 0.01. Hence, the intronic sequences are aligned just for this quality control step where sequences not passing the abovementioned criteria are filtered out and the sequences that pass are stored in files with extension '.filtered'. The custom scripts from Alam et al. were used for this step. It additionally yeilds a file where the summary of the changes made to each sequence file are logged. The output is -

```
FILE	REF	Nbefore	Nafter	total_distance	stats
B0348.4d.1_exon_0_83_DownstreamIntronFragment_84_154.fasta	CELEG.B0348.4d	34	50	165.0	too close:16 
B0348.4d.1_exon_10333_10575_DownstreamIntronFragment_10576_10646.fasta	CELEG.B0348.4d	36	50	175.0	too close:14 
B0348.4d.1_exon_14008_14149_DownstreamIntronFragment_14150_14220.fasta	CELEG.B0348.4d	40	50	195.0	Ns:1 too close:9 
B0348.4d.1_exon_14008_14149_UpstreamIntronFragment_13937_14007.fasta	CELEG.B0348.4d	43	46	210.0	too close:3 
B0348.4d.1_exon_1488_1849_DownstreamIntronFragment_1850_1920.fasta	CELEG.B0348.4d	45	50	220.0	too close:5 
B0348.4d.1_exon_1488_1849_UpstreamIntronFragment_1417_1487.fasta	CELEG.B0348.4d	42	46	205.0	short:1 too close:3 
```

Once the sequences are verified and filtered, we fix the orientation. For the downstream analysis, all intronic sequences, irrespective of the original strandedness of the gene, are converted to positive strand orientation. The shell script revcomp.sh is used which in turn calls a python script called revcomp_afterIntronFragments.py.finalintrons_aln is the name of the folder containing the alignments of the extracted intron fragments.

```
nohup bash revcomp_Introns.sh ./Positive_Strand_gene_names.txt + finalintrons_aln ./output_folder/ &
nohup bash revcomp_Introns.sh ./Negative_Strand_gene_names.txt - finalintrons_aln ./output_folder/ &
```

# References

Katoh, K., Rozewicki, J., & Yamada, K. D. (2019). MAFFT online service: multiple sequence alignment, interactive sequence choice and visualization. 20(4), 1160–1166. https://doi.org/10.1093/bib/bbx108

Alam, A., Duncan, A. G., Mitchell, J. A., & Moses, A. M. (biorxiv). Functional similarity of non-coding regions is revealed in phylogenetic average motif score representations. https://doi.org/10.1101/2023.04.09.536185


