from pymol import cmd

def clean(s = "all", renumber = True, ):
    selection_to_remove = s+ "and (organic or solvent)"
    cmd.remove(selection_to_remove)
