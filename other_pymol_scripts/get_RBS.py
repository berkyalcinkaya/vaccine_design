from pymol import cmd
from pymol_utils import remove_non_RBS
from utils import get_start_and_end
from io import StringIO
from Bio import SeqIO


def get_RBS(s, RBS_path = "data/consensus_RBS.fasta"):

    # MUST BE LOCATED IN VACCINE_DESIGN DIR TO USE
    RBS_seq_file = RBS_path
    RBS_record = SeqIO.parse(RBS_seq_file, "fasta")
    RBS_seq = (next(RBS_record)).seq   

    #isolate RBS for each object to align
    object_list = cmd.get_names(s)
    for object in object_list:
        RBS_object = object+"_RBS"
        cmd.copy(object, RBS_object )
        pymol_fasta_string = cmd.get_fastastr(RBS_object)
        fasta_io = StringIO(pymol_fasta_string) 
        sample_record = next(SeqIO.parse(fasta_io, "fasta"))
        sample_seq = sample_record.seq
        start,end = get_start_and_end(sample_seq, RBS_seq)
        remove_non_RBS(RBS_object,start,end)