'''
Berk Yalcinkaya
2/9/22

Pymol script for loading IRB structure files from tsv file
'''

import os
from glob import glob
import pandas as pd
from pymol import cmd


def load_all(fname):
    df = pd.read_csv(fname, sep='\t')
    pdb_ids = df["PDB ID"].tolist()
    for pdb_id in pdb_ids:
        cmd.fetch(pdb_id)

cmd.extend("load_all", load_all)
