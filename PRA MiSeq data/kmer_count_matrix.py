import os
import sys

def generate_dna_kmers(k):
    '''
    Return a list of all possible substrings of
    length k using only characters A, C, T, and G
    '''
    bases = ["A", "C", "U", "G"]

    last = bases
    current = []
    for i in range(k-1):
        for b in bases:
            for l in last:
                current.append(l+b)
        last = current
        current= []
    return last

input_PRA = sys.argv[1]
kmer_length = int(sys.argv[2])
output_kmer_count_matrix = sys.argv[3]
print("\n********************\n")
print("\n\nThe kmer length supplied by user : %d\n" %kmer_length)
print("Generating kmers of length %d from the PRA elements\n" %kmer_length)
print("\n********************\n")
kmers = generate_dna_kmers(kmer_length)

with open('All_kmers.txt', 'w') as handle:
    for f in kmers:
        handle.write(f)
        handle.write('\n')
handle.close()

PRALis=[]
with open(input_PRA, 'r') as r:
    for a in r:
        PRALis.append(str(a).strip())
print("Total number of PRA elements in the file : %d" %len(PRALis))

dic={}
i=0
import pandas as pd
for a in PRALis:
    Count = []
    for b in kmers:
        Count.append(a.count(b))
    dic.update({a:Count})
df = pd.DataFrame.from_dict(dic)

print("The %d-count matrix has been written in the output file! Exiting now" %kmer_length)
