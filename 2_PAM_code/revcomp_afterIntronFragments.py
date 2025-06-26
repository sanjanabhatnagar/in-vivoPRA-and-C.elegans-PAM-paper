import sys
from Bio import AlignIO
import os


input_alnfasta = sys.argv[1]
strand = sys.argv[2]
output_alnfasta = sys.argv[3]

alignment = AlignIO.read(input_alnfasta, "fasta")

def complement(s):
    l = []
    for a in s:
        if a == 'A':
            l.append('T')
        elif a == 'a':
            l.append('t')
        elif a == 'G':
            l.append('C')
        elif a == 'g':
            l.append('c')
        elif a == 'C':
            l.append('G')
        elif a == 'c':
            l.append('g')
        elif a == 'T':
            l.append('A')
        elif a == 't':
            l.append('a')
        else:
            l.append('')
    string = ''.join(l)
    return string

def reverse_complement(func, z):
    return func(z)[::-1]


type_aln = ''
for record in alignment:
    record_id = str(record.id)
    print(record_id)
    if (record_id.startswith('_R_CELEG')) and (strand == '+'):
        type_aln = 'strn_change'
        print(f"{record_id} in {strand} file {input_alnfasta}\n")
    elif (record_id.startswith('CELEG')) and (strand == '-'):
        type_aln = 'strn_change'
        print(f"{record_id} in {strand} file {input_alnfasta}\n")
    else:
        continue

if type_aln == 'strn_change':
    with open(output_alnfasta, 'a') as handle:
        for record in alignment:
            record_id = record.id
            seq_i = record.seq.reverse_complement()
            if '_R_' not in record.id:
                record_id_modified = f"_R_{record.id}"
            else:
                record_id_modified = record.id.replace('_R_', '')
            handle.write(f">{record_id_modified}\n{seq_i}\n")
    handle.close()
    os.remove(input_alnfasta)
