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
    '''removes non amino acid atoms and sets all chains to A if more than 1 chain is present'''
    if len(cmd.get_chains(s))>1:
        print("altering chains")
        cmd.alter(s, "chain = 'A'")
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
    print("consensus aln ", aln_consensus)
    print("sample aln ", aln_sample)
    return aln_consensus, aln_sample


def generate_edges(min, max, n = 8):
    bw =  (max - min)/n # bin width = bw
    edges = np.arange(min+bw, max+bw, bw)
    return edges


def get_color_name(entropy_val, edges):
    e = edges.tolist()
    e.reverse()
    i = 0
    while i<len(e) and entropy_val<=e[i]:
        i+=1
    return f"color_{i-1}"


def color_by_conservation(a, s = "all", subtypes = [1], save = False, thresh = 1.25):
    # define colors in pymol
    colors = [ [1,1,1], [0.63,0.63,0.63] , [0.50,0.68,0.56] , 
          [0.26,0.77,0.44] , [0.12,0.83,0.37] , 
          [0.01,0.89,0.31] , [0,0.95,0.30], [0.,1,0.70], [0,1,1] ]
    for i,color in enumerate(colors):
        color_string = f"color_{i}"
        cmd.set_color(color_string,color)

    # remove non amino acids
    clean_sequence(s)

    print(type(subtypes))
    
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
    cmd.iterate("(name ca)","stored.sample_resi_list.append(resi)") # append residue number
    print(len(stored.sample_resi_list))
    #print(stored.sample_resi_list)
    
    # produce global alignments
    consensus_aln, sample_algn = align(consensus_seq, sample_seq)
    
    # load consensus data
    # residue frequencies is stored in H{n}_H{m}_info.csv
    header = get_consensus_file_header(subtypes)
    consensus_info_path =  "data/" + header + "entropy.csv"
    entropy_data = pd.read_csv(consensus_info_path)

    # produces edges based on range of entropy values in entropy.csv
    # edges result in a number of bins equal to number of colors
    edges = generate_edges(entropy_data["entropy"].min(), 
                   entropy_data["entropy"].max(), 
                   n = len(colors)
                  )
    #print('edges\n:', edges)
    # must have same dimensions
    # colors are indentified by the upper bounding edge for entropy values
    assert len(edges) == len(colors)

    # produce the legend and save it
    legend = generate_legend(edges, colors, header[0:-1]) # last character includes an underscore,
                                                          # which we don't want in legend title
    legend_name = f"{header}color_legend.png"
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
            continue

        # perform coloring based on aa frequency at that position
        elif sample_algn[idx]!="-":
            sample_aa = get_sample_res(idx) # find aligned res in sample
            entropy = entropy_data.iloc[ref_aa_pos]["entropy"]
            color_name = get_color_name(entropy, edges)
            print("consenus aa, pos:", aa, ref_aa_pos+1, "|sample aa:", sample_aa, "|color:", color_name, "|entropy:", entropy)
            selection_string = f"resi {sample_aa}" 
            cmd.select("selection",selection_string)
            cmd.color(color_name, "selection")
            
            # increment reference sequence aa position only if aa is not a gap char
            ref_aa_pos+=1
        else:
            ref_aa_pos+=1

# add this function to pymol memory
cmd.extend("color_by_conservation", color_by_conservation)
