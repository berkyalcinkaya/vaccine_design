from pymol import cmd

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