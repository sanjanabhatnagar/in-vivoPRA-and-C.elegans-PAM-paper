# Following code was used for muscle and neuronal plasmid library for consistency
import pandas as pd
def complement(s):
    l = []
    for a in s:
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

f = open("~/PL1_merged.fastq", 'r').readlines() # Can be similarly run for all plasmid pools from 2-5
o1 = open("~/PL1_elements_and_flankingNt_region.txt", 'w') # Regions supplied to El_Extract() 
o3 = open("~/PL1_DiscardedSeqs.csv",'w') # Do not fit any of the flanking nts. matching criteria
o4 = open("~/PL1_Cis-Element_Barcodes_UMI.csv", 'w') # Extracted regions are processed using El_Extract()
o1.write('Index\tCis-element\tBarcode\tUMI\n')
o4.write('Cis-element,Barcode,UMI\n')

def El_Extract(x,y,z):
    h = (x[x.find('CAAATTGG')+len('CAAATTGG'):x.find('TGCTATGTCGTT')])
    k = (y[y.find('CGTTAA')+len('CGTTAA'):y.find('TACGTACTT')])
    i = (z[z.find('CTCTTC')+len('CTCTTC'):z.find('GACTT')])
    return [h,k,i]
i=0 
l=0
df = pd.DataFrame(columns=['Element_combinations'])
Array=[]
for a in f:
    a = a.translate(str.maketrans('', '', '[]"\' ,'))
    if a.startswith('CTAAACT'):
        cis = a[86:140]
        Barcode = a[-101:-47]
        UMI = a[-63:-11] 
        o1.write('%d\t%s\t%s\t%s\n' % (i, cis, Barcode, UMI))
        o4.write('%s\n' % (El_Extract(cis,Barcode,UMI)))
        i = i+1
    elif a.startswith('CAACAGT'):
        s = reverse_complement(complement, a)
        cis = s[86:140]
        Barcode = s[-101:-47]
        UMI = s[-63:-11]  
        o1.write('%d\t%s\t%s\t%s\n' % (i, cis, Barcode, UMI))
        o4.write('%s\n' % (El_Extract(cis,Barcode,UMI)))
        i = i+1
    else:
        l = l+1
        o3.write('%d\t%s' %(l, a))

o1.close()
o3.close()
