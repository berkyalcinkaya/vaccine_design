from io import StringIO
from os import remove
from Bio import SeqIO, pairwise2
import pandas as pd
import numpy as np
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import sys
import os

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
