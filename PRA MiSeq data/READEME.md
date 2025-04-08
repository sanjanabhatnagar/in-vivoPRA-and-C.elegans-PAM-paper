As described in the methods section, the quality scores of sequences generated from paired-end MiSeq multiplexed sequencing are assessed. The R1 and R2 files contain reads from both ends of the cDNA sequences. Before further analysis, these paired-end reads are merged into a single file to ensure both ends contribute to the calculation of PSI values for the corresponding reporter they map back to. This is followed by classifying the reads into different categories such as -
1. Spliced in isoforms
2. Spliced out isoforms
3. Unspliced isoforms
4. Cryptic spliced isoforms
5. Unknown isoforms

The script Splicing_Classifier_PRA_PSIcal.py is used to classify the reads as well as calculate the PSI values by collapsing reads using UMIs for each reporter, depicted by a unique barcode (which was introduced downstream of exon 3 in the minigene reporter and is thus present in the spliced isoforms) in the dataset. The script can be run as follows - 

```
python Splicing_Classifier_PRA_PSIcal.py CL1_merged.fastq
python Splicing_Classifier_PRA_PSIcal.py CL2_merged.fastq
```
