import sys
import os
from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
import pandas as pd

parent_name = sys.argv[1]
alignment_file = sys.argv[2] # Will take the alignment file with extension .aln 
reference_id = "CELEG"  # ID I have in the file

#Following files are generated from gff_exon-intron_annotation.py 
#github repository - https://github.com/sanjanabhatnagar/Inferring-Exon-and-Intron-Metadata-from-.gff-file.git
ExonIntron_coordinartesFile = sys.argv[3]
coordinate_df_exonintron = pd.read_csv(ExonIntron_coordinartesFile, sep='\t')

flank_size=None # Optional, I didn't use it for Cels since my gff_exon-intron_annotation.py script lets me select what size fragments I want. So my data is ready to use.
filtered_df = coordinate_df_exonintron[coordinate_df_exonintron['Parent'].str.contains(parent_name)][['Parent','group','ID','start','end','strand']]

output_dir = sys.argv[4]  # Specify the directory for intron fragments.

extract_gapped_exon_intron_sequences(alignment_file, reference_id, filtered_df, flank_size, output_dir)

def find_sps_fragment_in_sps_seq(sps_aln_fragment, sps_wholegene_aln):
    aln_fragmn_positionStart = sps_wholegene_aln.find(sps_aln_fragment)
    aln_fragmn_positionEnd = aln_fragmn_positionStart + len(sps_aln_fragment) + 1
    return(aln_fragmn_positionStart, aln_fragmn_positionEnd)

def map_ungapped_to_gapped(aligned_seq, ungapped_start, ungapped_end):
    gapped_start = None
    gapped_end = None
    ungapped_pos = 0

    for i, char in enumerate(aligned_seq, start=0):
        if char != '-':
            if ungapped_pos == ungapped_start: 
                gapped_start = i
            if ungapped_pos == ungapped_end:
                gapped_end = i
                break
            ungapped_pos += 1

    if gapped_end is None:
        gapped_end = len(aligned_seq)

    return gapped_start, gapped_end


def extract_gapped_exon_intron_sequences(alignment_file, reference_id, coords_df,flank_size, output_dir):
    coords_df.sort_values(by='start', ascending=True, inplace=True)
    coords_df['size']=abs(coords_df['start']-coords_df['end'])
    coords_exonsdf = coords_df[coords_df['group']=='exon']
    strand = coords_df['strand'].unique()[0]
    exon_coords = list(coords_exonsdf.apply(lambda row: (row['start'], row['end']), axis=1))

    alignment = AlignIO.read(alignment_file, "fasta")

    # Find the reference species sequence in the alignment
    reference_seq = None
    for record in alignment:
        if reference_id in record.id:
            reference_seq = record
            break
    if reference_seq is None:
        raise ValueError(f"Reference species {reference_id} not found in alignment.")

    seen_introns = set()
    for row in coords_exonsdf.itertuples(index=False):
        exon_start = row.start
        exon_end = row.end
        print(f'The exon start : {exon_start}. The exon end : {exon_end}')
        gapped_start, gapped_end = map_ungapped_to_gapped(reference_seq.seq, exon_start, exon_end)

        # If I want to do further resizing, I can add these lines of code.
        if flank_size==None:
            start = gapped_start
            end = gapped_end
        else:
            start = max(0, gapped_start - flank_size)
            end = min(len(reference_seq.seq), gapped_end + flank_size)
            # Collect SeqRecord fragments for the current exon region in C.elegans, Now in that seq record I need to
            # find the starts and ends for all other species.

        exon_alignment = MultipleSeqAlignment([SeqRecord(record.seq[start:end], id=record.id, description="")
            for record in alignment])


        # Get ungapped coordinates for the upstream intron from coords_df
        # So here, I'd have to work for CELEG and _R_CELEG
        if strand == '+':

            intron_coords = coords_df.loc[
                    ((coords_df['start'] == exon_end + 1) | (coords_df['end'] == exon_start - 1)) &
                    (coords_df['group'] == 'intron'),
                    ['start', 'end']
                    ]

        elif strand == '-':
            intron_coords = coords_df.loc[
                    ((coords_df['end'] == exon_start - 1) | (coords_df['start'] == exon_end + 1)) &
                    (coords_df['group'] == 'intron'),
                    ['start', 'end']
                    ]

        print(intron_coords)

        if intron_coords.empty:
            continue

        for _, intron_row in intron_coords.iterrows():
            CelsIntronStart, CelsIntronEnd = intron_row['start'], intron_row['end']
            cels_intron_length = abs(CelsIntronStart - CelsIntronEnd)
            print(f"The length of C.elegans intron fragment is : {cels_intron_length}")
            print(f'The relative ungapped coordinates are {CelsIntronStart} and {CelsIntronEnd}')

            intron_key = (CelsIntronStart, CelsIntronEnd)
            if intron_key in seen_introns:
                continue
        # Mapping the ungapped coordinates of introns to gapped coordinates, only for C.elegans!
            gapped_intron_start, gapped_intron_end = map_ungapped_to_gapped(reference_seq.seq, CelsIntronStart, CelsIntronEnd)
            print(f'The corresponding gapped coordinates are {gapped_intron_start} and {gapped_intron_end}')
            print('*-------------------------------------*-*-------------------------------------*')

            if (CelsIntronEnd == (exon_start - 1)):
                if strand == '+':
                    output_file = f"{output_dir}/{parent_name}_exon_{exon_start}_{exon_end}_UpstreamIntronFragment_{CelsIntronStart}_{CelsIntronEnd}.fasta"
                    frgmn_type="up"
                elif strand == '-':
                    output_file = f"{output_dir}/{parent_name}_exon_{exon_start}_{exon_end}_DownstreamIntronFragment_{CelsIntronStart}_{CelsIntronEnd}.fasta"
                    frgmn_type="down"

            elif (CelsIntronStart == (exon_end + 1)):
                if strand == '+':
                    frgmn_type="down"
                    output_file = f"{output_dir}/{parent_name}_exon_{exon_start}_{exon_end}_DownstreamIntronFragment_{CelsIntronStart}_{CelsIntronEnd}.fasta"
                elif strand == '-':
                    frgmn_type="up"
                    output_file = f"{output_dir}/{parent_name}_exon_{exon_start}_{exon_end}_UpstreamIntronFragment_{CelsIntronStart}_{CelsIntronEnd}.fasta"

            with open(output_file, "w") as intron_handle:
                Cels_intron_ungapped = ""
                record_intron = ""
                reference_id_aln = None
                record_intron_Cels = ""
          
                for record in alignment:
                    record_id = record.id
                    j = 0
                    for aln_fragment in exon_alignment:
                        if reference_id in record_id:
                            # taking intron based on the gff inferred coordinates and size from Celegans aligned seq
                            reference_seq_aln = record.seq
                            reference_id_aln = record.id
                            Cels_intron_gapped = reference_seq_aln[gapped_intron_start:gapped_intron_end] # Here it takes both upstream and downstream intron correctly for C.elegans/ref species
                            while (len(Cels_intron_ungapped) <= cels_intron_length) and (j < len(Cels_intron_gapped)):
                                if Cels_intron_gapped[j] != '-':
                                    Cels_intron_ungapped += Cels_intron_gapped[j]
                                j += 1
                            record_intron_Cels = Cels_intron_ungapped #And here I have correct sequence for a given intron upstream or downstream.

                        elif record_id == aln_fragment.id:
                        # Locate aln_fragment in record.seq and get ungapped start and ungapped end.
                            sps_start, sps_end = find_sps_fragment_in_sps_seq(aln_fragment.seq, record.seq)
                            intron_frm_aln = ''  # To store the extracted sequence
                            count = 0  # Counter for non-gap nucleotides

                            if strand == '-':
                                if frgmn_type == "down":
                                    i = sps_start - 1
                                    while count < cels_intron_length and i >= 0:
                                        if record.seq[i] != '-':
                                            intron_frm_aln += record.seq[i]
                                            count += 1
                                        i -= 1

                                if frgmn_type == "up":
                                    i = sps_end + 1
                                    while count < cels_intron_length and i < len(record.seq):
                                        if record.seq[i] != '-':
                                            intron_frm_aln += record.seq[i]
                                            count += 1
                                        i += 1

                            elif strand == '+':
                                if frgmn_type == "down":
                                    i = sps_end + 1
                                    while count < cels_intron_length and i < len(record.seq):
                                        if record.seq[i] != '-':
                                            intron_frm_aln += record.seq[i]
                                            count += 1
                                        i += 1

                                elif frgmn_type == "up":
                                    i = sps_start - 1
                                    while count < cels_intron_length and i >= 0:
                                        if record.seq[i] != '-':
                                            intron_frm_aln += record.seq[i]
                                            count += 1
                                        i -= 1

                            record_intron = intron_frm_aln
                            intron_handle.write(f">{record.id}\n{record_intron}\n")

                intron_handle.write(f">{reference_id_aln}\n{record_intron_Cels}\n")

~                                                                                                               
