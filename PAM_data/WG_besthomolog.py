import sys
from Bio import SeqIO

def best_homolog_wholegene(input_aa_file, input_fa_file, output_fa_file):
    multiple_aa_aln = SeqIO.parse(input_aa_file, "fasta")
    seq_file_list = list(SeqIO.parse(input_fa_file, "fasta"))

    with open(output_fa_file, "w") as out_f:
        for aa in multiple_aa_aln:
            aa_id = aa.id.split('_PID')[0]
            print(f"The amino acid id {aa_id} is being searched")
            SeqIO.write([seq for seq in seq_file_list if aa_id in seq.id], out_f, "fasta")

input_file = sys.argv[1]
input_2_file = sys.argv[2]
output_file = sys.argv[3]

best_homolog_wholegene(input_file, input_2_file, output_file)
