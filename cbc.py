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
import numpy as np
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from io import StringIO
from os import remove
from Bio import SeqIO, pairwise2
import pandas as pd
import numpy as np
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import sys
import os


def clean_sequence(s):
    '''removes non amino acid atoms and sets all chains to A if more than 1 chain is present'''
    if len(cmd.get_chains(s))>1:
        print("altering chains")
        cmd.alter(s, "chain = 'A'")
    selection_to_remove = s+ " and (organic or solvent)"
    cmd.remove(selection_to_remove)

def remove_non_RBS(s, start, end):
    # pymol selection is 1 indexed and inclusive
    to_remove = f"{s} and not resi {start+1}-{end+1}"
    cmd.remove(to_remove)


def local_align(consensus, RBS, verbose = True):
    '''peforms a global alignment between a consensus and sample sequence
    gap open penalty is raised very high to prevent gaps in consensus sequence'''
    # def scoring metrics
    gap_open_default = -3
    gap_open_consensus = -3
    gap_extend_default = -0.1
    alns = pairwise2.align.localxd(sequenceA=consensus, sequenceB=RBS, 
                                      openA = gap_open_consensus, 
                                      extendA = gap_extend_default, 
                                      openB = gap_open_default, 
                                      extendB = gap_extend_default)
    top_aln = alns[0] # choose first (best) alignment
    aln_consensus, aln_RBS, _, _, _ = top_aln
    
    if verbose:
        print("consensus aln ", aln_consensus)
        print("RBS aln ", aln_RBS)
    
    return aln_consensus, aln_RBS


def get_start_and_end(consensus, RBS, verbose = True):
    '''gets the beginning and end residue of the RBS in a given HA sequence'''
    aln_consensus, aln_RBS = local_align(consensus, RBS, verbose = verbose)
    
    start = -1
    first_char = False
    for aa_char in aln_RBS:
        start +=1
        if aa_char!="-":
            break
    
    end = len(aln_RBS)
    for aa_char in reversed(aln_RBS):
        end -=1
        if aa_char!="-":
            break
    
    if verbose:
        print("Alignment Indexes ", "start: ", start, "end: ", end)
        print(aln_consensus[start:end+1])
    
    return aln_idx_to_consensus_idx(start, end, aln_consensus)


def aln_idx_to_consensus_idx(start, end, consensus_aln):
    '''helper method for get_start_and_end. Converts the alignment index to a normal sequence index'''
    idx = -1
    for char in consensus_aln[0:start+1]:
        if char!="-":
            idx+=1
    start_idx = idx
    
    for char in consensus_aln[start+1:end+1]:
        if char!="-":
            idx+=1
    end_idx = idx
    
    return start_idx,end_idx 


def color_by_conservation(a, s = "all", subtypes = "1", save = 1):
    '''HOW TO USE
    Because of the way pymol accepts arguments, you must always pass the argument a as a dummy argument when you 
    run the program. In other words, you must call 
    color_by_conservation <any integer> 
    to use this program with default arguments
    
    Examples:
    H1 generation
    --------------
    color_by_conservation 1, subtypes = 1

    H1-H3 generation
    --------------
    color_by_conservation 1, subtypes = 1,3

    H2 generation
    --------------
    color_by_conservation 1, subtypes = 2
    
    '''
    if save>=1:
        save = True
    elif save == 0:
        save = False

    # define colors in pymol
    colors = [ [1,1,1], [0.63,0.63,0.63] , [0.50,0.68,0.56] , 
          [0.26,0.77,0.44] , [0.12,0.83,0.37] , 
          [0.01,0.89,0.31] , [0,0.95,0.30], [0.,1,0.70], [0,1,1] ]
    for i,color in enumerate(colors):
        color_string = f"color_{i}"
        cmd.set_color(color_string,color)

    # remove non amino acids
    clean_sequence(s)
    
    # load consensus sequence using Biopython
    # most consensus sequences can be accessed by the subtype number 
    # ex: 1 in H1 or 1,3 in H1_H3
    if subtypes == "human":
        file_header = "human_"
    else:
        file_header = get_consensus_file_header(subtypes)
    consensus_seq_file = "data/" + file_header + "consensus.fasta"
    consensus_record = SeqIO.parse(consensus_seq_file, "fasta")
    consensus_seq = next(consensus_record).seq

    # get loaded pymol structure's sequence
    # use Biopython and StringIO to load the string directly to a SeqIO object
    pymol_fasta_string = cmd.get_fastastr(s)
    fasta_io = StringIO(pymol_fasta_string) 
    sample_record = next(SeqIO.parse(fasta_io, "fasta"))
    sample_seq = sample_record.seq
    print(sample_seq)

    # populate residue list
    # name of all residues are now accessible and can be accessed by integer index
    # irregular numbering naming is handled in this way
    stored.sample_resi_list=[]
    cmd.iterate(f"({s} and name ca)","stored.sample_resi_list.append(resi)") # append residue number
    #print(len(stored.sample_resi_list))
    print(stored.sample_resi_list)
    
    # produce global alignments
    consensus_aln, sample_algn = align(consensus_seq, sample_seq)
    
    # load consensus data
    # residue frequencies is stored in H{n}_H{m}_info.csv
    consensus_info_path =  "data/" + file_header + "entropy.csv"
    entropy_data = pd.read_csv(consensus_info_path)

    # produces edges based on range of entropy values in entropy.csv
    # edges result in a number of bins equal to number of colors
    edges = generate_edges(entropy_data["entropy"].min(), 
                   entropy_data["entropy"].max(), 
                   n = len(colors)
                  )

    # must have same dimensions
    # colors are indentified by the upper bounding edge for entropy values
    print("Min Entropy in HA: ", entropy_data["entropy"].min(), "Max Entropy: ", entropy_data["entropy"].max())
    #print(edges)
    assert len(edges) == len(colors)

    # produce the legend and save it
    legend = generate_legend(edges, colors, file_header[0:-1]) # last character includes an underscore,
                                                          # which we don't want in legend title
    if save:
        legend_name = f"{file_header}color_legend.png"
        print("\nsaving entropy legend to ",legend_name)
        export_legend(legend, filename = legend_name)


    ## Helper Function for Mapping Consensus onto Experimental Structure
    #-----------------------------------------------------------------------
    def get_sample_res(ref_idx):
        '''given an index of an aa from the aligned reference sequence, 
        find the pymol resi number that can be used to select the amino acid'''
        i=-1
        for char in sample_algn[0:ref_idx+1]:
            if char!="-":
                i+=1
        pymol_resi = stored.sample_resi_list[i] # get given residue name from resi_list
        return pymol_resi
    #----------------------------------------------------------------------


    ref_aa_pos = 0 # because of gaps in ref seq, idx is not the same as the aa pos
                    # keep track of aa pos with this variable
                    # Because of global alignment parameters, there SHOULDNT be any gaps
                    # Albeit, this counter exists as sanity check
    for idx,aa in enumerate(consensus_aln):
        # skip over gaps in reference sequence
        if aa=="-":
            print("gap in consensus")
            continue

        # perform coloring based on aa frequency at that position
        elif sample_algn[idx]!="-":
            sample_aa = get_sample_res(idx) # find aligned res in sample
            entropy = entropy_data.iloc[ref_aa_pos]["entropy"]
            color_name = get_color_name(entropy, edges)
            print("consenus aa, pos:", aa, ref_aa_pos+1, "|sample aa:", sample_aa, "|color:", color_name, "|entropy:", entropy)
            selection_string = f"{s} and resi {sample_aa}" 
            cmd.select("selection",selection_string)
            cmd.color(color_name, "selection")
            
            # increment reference sequence aa position only if aa is not a gap char
            ref_aa_pos+=1
        else:
            print("gap in sample")
            ref_aa_pos+=1
    #print(consensus_aln)



def color_by_conservation_rbs(subtypes, s = "all", save = True):
    '''Given an RBS predicted structure, colors residues based on entropy at that point and generates a figure
    
    USAGE
    color_by_conservation_rbs'''
    
    colors = [ [1,1,1], [0.63,0.63,0.63] , [0.50,0.68,0.56] , 
          [0.26,0.77,0.44] , [0.12,0.83,0.37] , 
          [0.01,0.89,0.31] , [0,0.95,0.30], [0.,1,0.70], [0,1,1] ]
    for i,color in enumerate(colors):
        color_string = f"color_{i}"
        cmd.set_color(color_string,color)

    # load consensus sequence using Biopython
    # most consensus sequences can be accessed by the subtype number 
    # ex: 1 in H1 or 1,3 in H1_H3
    if subtypes == "human":
        file_header = "human_"
    else:
        file_header = get_consensus_file_header(subtypes)
    consensus_seq_file = "data/" + file_header + "consensus.fasta"
    consensus_record = SeqIO.parse(consensus_seq_file, "fasta")
    consensus_seq = next(consensus_record).seq  
    
    # get loaded pymol structure's sequence
    # use Biopython and StringIO to load the string directly to a SeqIO object
    pymol_fasta_string = cmd.get_fastastr(s)
    fasta_io = StringIO(pymol_fasta_string) 
    RBS_record = next(SeqIO.parse(fasta_io, "fasta"))
    RBS_seq = RBS_record.seq

    RBS_start, RBS_end = get_start_and_end(consensus_seq, RBS_seq, verbose = True)
    print("RBS start: ", RBS_start, "| RBS End: ", RBS_end)

    # populate residue list
    # name of all residues are now accessible and can be accessed by integer index
    # irregular numbering naming is handled in this way
    stored.sample_resi_list=[]
    cmd.iterate(f"({s} and name ca)","stored.sample_resi_list.append(resi)") # append residue number
    #print(len(stored.sample_resi_list))
    #print(stored.sample_resi_list)

    # load consensus data
    # residue frequencies is stored in H{n}_H{m}_info.csv
    consensus_info_path =  "data/" + file_header + "entropy.csv"
    entropy_data = pd.read_csv(consensus_info_path).iloc[RBS_start:RBS_end+1]
    #print(entropy_data)

      # produces edges based on range of entropy values in entropy.csv
    # edges result in a number of bins equal to number of colors
    edges = generate_edges(entropy_data["entropy"].min(), 
                   entropy_data["entropy"].max(), 
                   n = len(colors)
                  )

    # must have same dimensions
    # colors are indentified by the upper bounding edge for entropy values
    print("Min Entropy: ", entropy_data["entropy"].min(), "Max Entropy: ", entropy_data["entropy"].max())
    #print(edges)
    assert len(edges) == len(colors)

    # produce the legend and save it
    legend = generate_legend(edges, colors, file_header[0:-1]) # last character includes an underscore,
                                                          # which we don't want in legend title
    if save:
        legend_name = f"{file_header}color_legend_RBS.png"
        print("\nsaving entropy legend to ",legend_name)
        export_legend(legend, filename = legend_name)

    for idx,aa in enumerate(RBS_seq):
        pymol_resi = stored.sample_resi_list[idx]
        entropy = entropy_data.iloc[idx]["entropy"]
        color_name = get_color_name(entropy, edges)
        print("RBS pos:", idx+1, "|RBS aa:", aa, "|color:", color_name, "|entropy:", entropy)
        selection_string = f"{s} and resi {pymol_resi}" 
        cmd.select("selection",selection_string)
        cmd.color(color_name, "selection")


############### Helper functions ##################

def generate_legend(edges, colors, legend_title):
    legend_handles = []
    i = 0
    for color, entropy_val in zip(reversed(colors), edges):
        lower_bound = 0 if i==0 else edges[i-1]
        s = mlines.Line2D([], [], color=tuple(color), marker='s', 
                        linestyle='None',linewidth = 20, markersize= 14, 
                        label=f"{round(lower_bound,2)} ≤ entropy < {round(entropy_val,2)}" )
        legend_handles.append(s)
        i+=1
    legend = plt.legend(handles = legend_handles, title = legend_title, 
                        title_fontsize = 15, facecolor = "black", 
                        labelcolor = "white", mode = "expand")
    plt.setp(legend.get_title(), color='white')
    
    # show legend 
    plt.axis("off")
    plt.show(block = False) # allow program to proceed with block=false

    return legend


def export_legend(legend, filename="legend.png"):
    fig  = legend.figure
    fig.canvas.draw()
    bbox  = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(filename, dpi="figure", bbox_inches=bbox)


def get_consensus_file_header(subtypes: str):
    '''pattern
    -----------
    subtypes = "1,2,3" --> file_header = HX_HY_HZ_'''
    
    if len(subtypes) == 1:
        return f"H{subtypes}_"
        
    else: 
        subtypes = subtypes.split(",")
        for subtype in subtypes.sort():
            file_string+=f"H{subtype}_"
        return file_string


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
    print("consensus aln ", aln_consensus)
    print("sample aln ", aln_sample)
    return aln_consensus, aln_sample


def generate_edges(min, max, n = 8):
    bw =  (max - min)/n # bin width = bw
    edges = np.linspace(min, max, n+1)
    return edges[1:]


def get_color_name(entropy_val, edges):
    e = edges.tolist()
    e.reverse()
    i = 0
    while i<len(e) and entropy_val<=e[i]:
        i+=1
    return f"color_{i-1}"

# add this function to pymol memory
cmd.extend("color_by_conservation", color_by_conservation)
cmd.extend("color_by_conservation_rbs", color_by_conservation_rbs)