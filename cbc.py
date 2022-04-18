'''
Berk Yalcinkaya
Given a user-inputted subtype, this pymol script accesses consensus sequence data from the 
/data folder and "maps" consensus sequence onto current pymol objects, coloring
the sample residues based on sequence entropy on a gradient from dark green (low entropy, high conservation)
to white (high entropy, above a threshold of 1.25). '''


from pymol import cmd,stored
from io import StringIO
from Bio import SeqIO, pairwise2
import pandas as pd


def get_consensus_file_header(subtypes:list):
    '''pattern
    -----------
    subtypes = [HX, HY, HZ] --> file_header = HX_HY_HZ_'''
    file_string = ""
    for subtype in reversed(subtypes):
        file_string+=f"H{subtype}_"
    #file_string+= "consensus.fasta"
    return file_string


def clean_sequence(s):
    '''removes '''
    selection_to_remove = s+ " and (organic or solvent)"
    cmd.remove(selection_to_remove)


def align(consensus, sample):
    '''peforms a global alignment between a consensus and sample sequence
    gap open penalty is raised very high to prevent gaps in consensus sequence'''
    # def scoring metrics
    gap_open_default = -10
    gap_open_consensus = -50
    gap_extend_default = -0.5
    alns = pairwise2.align.globalxd(sequenceA=consensus, sequenceB=sample, 
                                      openA = gap_open_consensus, 
                                      extendA = gap_extend_default, 
                                      openB = gap_open_default, 
                                      extendB = gap_extend_default)
    top_aln = alns[0] # choose first (best) alignment
    aln_consensus, aln_sample, _, _, _ = top_aln
    return aln_consensus, aln_sample


def color_by_conservation(s = "all", subtypes = [1], save = False, thresh = 1.25):

    # remove non amino acids
    clean_sequence(s)
    
    # load consensus sequence using Biopython
    consensus_seq_file = "data/" + get_consensus_file_header(subtypes)+"consensus.fasta"
    consensus_record = SeqIO.parse(consensus_seq_file, "fasta")
    consensus_seq = next(consensus_record).seq

    # get loaded pymol structure's sequence
    # use Biopython and StringIO to load the string directly to a SeqIO object
    pymol_fasta_string = cmd.get_fastastr("all")
    fasta_io = StringIO(pymol_fasta_string) 
    sample_record = next(SeqIO.parse(fasta_io, "fasta"))
    sample_seq = sample_record.seq

    # populate residue list
    # name of all residues are now accessible and can be accessed by integer index
    # irregular numbering naming is handled in this way
    stored.sample_resi_list=[]
    cmd.iterate("(name ca)","stored.list.append(resi)") # append residue number

    # check for irregular numbering and alert user
    for resi in stored.sample_resi_list:
        try: 
            int(resi)
        except:
            print("irregular residue numbering detected", resi)

    # produce alignments
    consensus_aln, sample_algn = align(consensus_seq, sample_seq)
    
    # load consensus data
    # residue frequencies is stored in H{n}_H{m}_info.csv
    consensus_info_path =  "data/" + get_consensus_file_header(subtypes) + "entropy.csv"
    consensus_data = pd.read_csv(consensus_info_path)

    def get_sample_res(ref_idx):
        '''given an index of an aa from the aligned reference sequence, 
        find the pymol resi number that can be used to select the amino acid'''
        idx=-1
        for char in sample_algn[0:ref_idx+1]:
            if char!="-":
                idx+=1
        pymol_resi = stored.sample_resi_list[idx] # get given residue name ffrom resi_list
        return pymol_resi

    ref_aa_pos = 0 # because of gaps in ref seq, idx is not the same as the aa pos
                    # keep track of aa pos with this variable
    for idx,aa in enumerate(consensus_aln):
        # skip over gaps in reference sequence
        if aa=="-":
            continue
        
        # perform coloring based on aa frequency at that position
        else:
            sample_aa = get_sample_res(idx)
            entropy = consensus_data.iloc[ref_aa_pos]["entropy"]
            
            # determine color to be used
            if entropy > thresh:
                color_name = 'white'
            else:
                entropy_score = entropy/thresh
                start_green = 550 # lighter green
                end_green = 615 # darker green, end of color spectrum
                # produce a color for the given resiude
                # lower entropy = lower entropy score = higher color number
                color = str(int((end_green - start_green)*(1-entropy_score) + start_green))
                color_name = 'c'+color # pymol one number coloring begins with c
        
            # select resiude and color it
            selection_string = f"resi {sample_aa}" 
            cmd.select("selection",selection_string)
            cmd.color(color_name, "selection")
            
            # increment reference sequence aa position only if aa is not a gap char
            ref_aa_pos+=1
# add this function to pymol memory
cmd.extend("color_by_conservation", color_by_conservation)
