from Bio import SeqIO
import sys

def calculate_pid(ref_seq, seq_frm_aln): # Moses lab function was used
    nid=0
    ndiff=0
    for pos in range(0,len(ref_seq)):
        if ('-' in ref_seq[pos]+seq_frm_aln[pos]):
            pass
        elif (ref_seq[pos] == seq_frm_aln[pos]):
            nid+=1 # number of identical pos (nts)
        else:
            ndiff+=1 # number of different nts at the same position

    if (nid == 0):
        return float(10)
    elif (ndiff == 0):
        return float(0)
    else:
        return float(ndiff)/float(nid+ndiff)


def select_best_ortholog(input_file, output_file, reference_id="CELEG"): # Reference species name can be changed here if running for a different organism
    sequences = list(SeqIO.parse(input_file, "fasta"))
    reference_seq = None
    best_matches = {}

    for seq in sequences:
        if reference_id in seq.id:
            reference_seq = seq.seq
            break

    if not reference_seq:
        raise ValueError(f"Reference {reference_id} not found in {input_file}")

    with open('PIDs_allparalogs.txt','a') as handle: # This file contains the IDs of all the orthologs from a given species and their PID with C.elegans ortholog in that ortholog group.
        handle.write(f'Event\tRef\tSps\tPID\n')
        for seq in sequences:
            species = seq.id.split(".")[0]
            if species == reference_id:
                continue

            pid = calculate_pid(reference_seq, seq.seq)
            handle.write(f'{input_file}\t{reference_id}\t{seq.id}\t{pid}\n')
            if species not in best_matches or pid > best_matches[species][1]:
                best_matches[species] = (seq, pid)  # Paralog with highest PID would get stored.

    with open(output_file, "w") as out_f:
        SeqIO.write([seq for seq in sequences if reference_id in seq.id], out_f, "fasta")
        for species, (seq,pid) in best_matches.items():

            seq.id = f"{seq.id}_PID:{pid}"
            seq.description = ""
            SeqIO.write(seq, out_f, "fasta")


input_file = sys.argv[1]
output_file = sys.argv[2]

select_best_ortholog(input_file, output_file)
