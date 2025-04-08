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
