'''
Berk Yalcinkaya
4/22/22
Pymol script that allows for viewing output from B and T cell epitope predictions on an structure
'''
def is_continuous(e):
    pass

def get_selection_string(e):
    pass

def visualize(epitopes: str, s = "all"):
    '''
    Inputs
    ---------------
    epitopes (str): a single quote (') sperated string of epitopes to be visualized in pymol
                    can be both linear and discontinuous
    
    example of acceptable discontinuous epitope: 
    "A:S100, A:D101, A:G102, A:N103, A:N106"

    Because we are using predicted structures, we will assume there is only a chain A

    Outputs
    --------------
    creates numbered selections for each epitopes with either _d or _c to indicate discontinous 
    or continouous, respectively
    '''
    # turn epitopes into a list
    epitopes = epitopes.split("'")
    epitopes =  [epitope.strip() for epitope in epitopes]

    for epitope in epitopes:
        pass

