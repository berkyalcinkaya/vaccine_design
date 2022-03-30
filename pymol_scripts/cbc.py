from pymol import cmd,stored
from io import StringIO
from Bio import SeqIO

# load consensus sequence using Biopython
consensus_seq_file = "H1_consensus.fasta"
consensus_record = SeqIO.parse(consensus_seq_file, "fasta")
consensus_seq = next(consensus_record).seq

# get loaded pymol structure's sequence
# use Biopython and StringIO to load the string directly to a SeqIO object
pymol_fasta_string = cmd.get_fastastr("all")
fasta_io = StringIO(pymol_fasta_string) 
sample_record = next(SeqIO.parse(fasta_io, "fasta"))
sample_seq = sample_record.seq

#populate residue list
stored.sample_resi_list=[]
cmd.iterate("(name ca)","stored.list.append(resi)") # append residue number

# check for irregular numbering
for resi in stored.sample_resi_list:
    try: 
        int(resi)
    except:
        print("irregular residue numbering detected", resi)



def get_sample_res(ref_idx):
    '''given an index of an aa from the aligned reference sequence, 
    find the pymol resi number that can be used to select the amino acid'''
    
    idx=-1
    for char in sample_seq[0:ref_idx+1]:
        if char!="-":
            idx+=1
    
    pymol_resi = stored.sample_resi_list[idx]
    return pymol_resi