import os
os.getcwd()
import sys 
import pandas as pd

#Functions
def complement(s):
    l = []
    for a in s:
        if a == ',':
            l.append(',')
        if a == 'A':
            l.append('T')
        elif a == 'G':
             l.append('C')
        elif a == 'C':
             l.append('G')
        elif a == 'T':
             l.append('A')
        else:
            l.append('')
    string = ''.join(l)
    return string

def reverse_complement(func, z):
    return func(z)[::-1]
    
def Hamming_check_0_or_1(Bar_to_map, bar_for_map):
    errors = 0
    x = len(bar_for_map)
    y = len(Bar_to_map)
    if x == y:
        for i in range(len(Bar_to_map)):
            if Bar_to_map[i] != bar_for_map[i]:
                errors += 1
        return errors
    elif y>x:
        l = 0
        m = x
        errors_lis=[]
        while (l<= (y-x)) and (m <= y):
            k_mer = Bar_to_map[l:m]
            l = l+1
            m = m+1
            errors_kmer=0
            for i in range(len(k_mer)):
                if k_mer[i] != bar_for_map[i]:
                    errors_kmer += 1
                i=i+1
            errors_lis.append(errors_kmer)
        return(min(errors_lis))
    else:
        return(4)

# CACATTCGTTAA - Exon3-Upstream region of barcode
# TACGTACTTCTGAG - Downstream region of barcode
# TACGTCTTCTGAG - Downstream region of barcode - has a deletion

def Splicing_Classifier(u,y,z):
# TTAAAGCTAA - Exon1-CE junction
# TGAAAAAAGAAGATCCA - CE-Exon3 junction
    if ((Hamming_check_0_or_1(u,'TTAAAGCTAA')<2)and(Hamming_check_0_or_1(u,'TGAAAAAAGAAGATCCA')<2)and(Hamming_check_0_or_1(u,'CACATTCGTTAA')<2)and((Hamming_check_0_or_1(u,'TACGTACTTCTGAG')<2)or(Hamming_check_0_or_1(u,'TACGTCTTCTGAG')<2))):
        k = (y[y.find('CGTTAA')+len('CGTTAA'):y.find('TACGT')])
        i = (z[z.find('ATTACTCTTC')+len('ATTACTCTTC'):z.find('CTTGCT')])
        return ('Included',k,i)
# TAAAGAAGAT - Exon1-Exon3 junction
# TGAAAAAAGAAGATCCA - Overkill, CE-exon3 junction should be absent
    elif ((Hamming_check_0_or_1(u,'TAAAGAAGAT')<2)and(Hamming_check_0_or_1(u,'CACATTCGTTAA')<2)and((Hamming_check_0_or_1(u,'TACGTACTTCTGAG')<2)or(Hamming_check_0_or_1(u,'TACGTCTTCTGAG')<2))and not(Hamming_check_0_or_1(u,'TGAAAAAAGAAGATCCA')<2)):
        k = (y[y.find('CGTTAA')+len('CGTTAA'):y.find('TACGT')])
        i = (z[z.find('ATTACTCTTC')+len('ATTACTCTTC'):z.find('CTTGCT')])
        return ('Skipped',k,i)
# TTAAAGCTAA - Exon1-CE junction
# GTTTCAAAT - Region between CE and introduced PRA element, region of intron2, 3nt upstream of PRA element
# AAATTGG  - If the read doesn't contain entire GTTTCAAAAT sequence, it may look for smaller AAATTGG region, flanks PRA element
# AAGATCCAT - should contain Exon 3
# GAAAAAAGAAGATCC - should not contain this region CE-Exon 3 junction
# TTTCAGAAGATC - should not contain this region Intron2-Exon3 junction
    elif ((Hamming_check_0_or_1(u,'TTAAAGCTAA')<2)and(Hamming_check_0_or_1(u,'AAGATCCAT')<2)and((Hamming_check_0_or_1(u,'GTTTCAAAT')<2)or(Hamming_check_0_or_1(u,'AAATTGG')<2))and(Hamming_check_0_or_1(u,'CACATTCGTTAA')<2)and((Hamming_check_0_or_1(u,'TACGTACTTCTGAG')<2)or(Hamming_check_0_or_1(u,'TACGTCTTCTGAG')<2))and not(Hamming_check_0_or_1(u,'GAAAAAAGAAGATCC')<2)and not(Hamming_check_0_or_1(u,'TTTCAGAAGATC')<2)): #'AAAAGGTGATCC not coming in reads, so checking for AAGATCCAT exon3 beginning in reads,
        k = (y[y.find('CGTTAA')+len('CGTTAA'):y.find('TACGT')])
        i = (z[z.find('ATTACTCTTC')+len('ATTACTCTTC'):z.find('CTTGCT')])
        return ('Cryptic',k,i)
# AATTTTTTAG - Either has some portion of intron 1
# CAAATTTTTTCAG - or has some portion of intron 2
    elif ((Hamming_check_0_or_1(u,'AATTTTTTAG')<2)or(Hamming_check_0_or_1(u,'CAAATTTTTTCAG')<2)and(Hamming_check_0_or_1(u,'CACATTCGTTAA')<2)and((Hamming_check_0_or_1(u,'TACGTACTTCTGAG')<2)or(Hamming_check_0_or_1(u,'TACGTCTTCTGAG')<2))):
        k = (y[y.find('CGTTAA')+len('CGTTAA'):y.find('TACGT')])
        i = (z[z.find('ATTACTCTTC')+len('ATTACTCTTC'):z.find('CTTGCT')])
        return ('Unspliced',k,i)
# AAAAGGTGATCC - CE-intron 2 junction
# TTTCAGAAGATC - intron2 - exon3 junction
    elif ((Hamming_check_0_or_1(u,'AAAAGGTGATCC')<2)and(Hamming_check_0_or_1(u,'TTTCAGAAGATC')<2)and(Hamming_check_0_or_1(u,'CACATTCGTTAA')<2)and((Hamming_check_0_or_1(u,'TACGTACTTCTGAG')<2)or(Hamming_check_0_or_1(u,'TACGTCTTCTGAG')<2))):
        k = (y[y.find('CGTTAA')+len('CGTTAA'):y.find('TACGT')])
        i = (z[z.find('ATTACTCTTC')+len('ATTACTCTTC'):z.find('CTTGCT')])
        return ('Unspliced',k,i)
    else:
        k = (y[y.find('CGTTAA')+len('CGTTAA'):y.find('TACGT')])
        i = (z[z.find('ATTACTCTTC')+len('ATTACTCTTC'):z.find('CTTGCT')])
        return ('Unknown',k,i)

    


#Extracting spliced-in isoforms
infile = sys.argv[1]
f_name, f_ext = os.path.splitext(Infile)
outfile1 = ("%s_SplicingInfo.tsv" %f_name)
outfile2 = ("%s_SplicedIsoformCounts.tsv" %f_name)

f = open(infile, 'r').readlines()
o4 = open(outfile1, 'w')


o4.write('Index\tCategory\tBarcode\tUMI\n')       
i=0 
l=0
Array=[]
for b in f:
    a = b.strip('[\"[\'')
    if (Hamming_check_0_or_1(a[:15],'GATTACAAGG')<3):
        UMI = a[-50:-10]
        Barcode = a[-88:-58]
        o4.write('%d\t%s\n' % (i, Splicing_Classifier(a, Barcode, UMI)))
        i = i+1
    elif (Hamming_check_0_or_1(a[:15],'CAACAGTACG')<3):
        UMI = reverse_complement(complement, (a[10:50]))
        Barcode = reverse_complement(complement,(a[48:78]))
        o4.write('%d\t%s\n' % (i, Splicing_Classifier(a, Barcode, UMI)))
        i = i+1
    else:
        l = l+1
        
o4.close()

r_df = pd.read_csv(outfile1,sep='\t')

df = pd.DataFrame(r_df)
df.drop(columns='Index', axis=0, inplace=True)

df['Barcode_Activity'] = (df['Barcode'] +'\t'+ df['Category']).astype(str)
df=df[['Barcode_Activity', 'UMI']]

df_remdup = pd.DataFrame.drop_duplicates(df)
(df_remdup['Barcode_Activity'].value_counts()).to_csv(outfile2)




