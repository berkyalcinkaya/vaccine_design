'''
4/15/22
Berk Yalcinkaya
This script runs the entire consensus sequence generation pipeline:
1) automated sequency query from IRD
2) cleaning, clustering, and aligning
3) consensus sequence generatior and entropy calculations using the protein-consensus-sequence 
    package

**IMPORTANT**
Add the path your muscle and cd-hit directory below. 
'''

import os
import subprocess
import pandas as pd
import shutil
import requests
from Bio import SeqIO
from io import StringIO
import sys

# GLOBALS
# ADD YOUR OWN PATHS HERE
#-----------------------------------------------
cd_hit_path = "/Users/berk/venvs/vaccine/bin/cd-hit" # change this
muscle_path = "/Users/berk/opt/muscle" # change this
#-----------------------------------------------
log_file = "log.txt"

def log(texts:list, newline = True, endline = True):
    '''write text to log_file path, assuming file has been created'''
    f = open(log_file, "a") # log_file is a global variable
    if newline:
        f.write("\n")
        
    for text in texts:
        f.write(text)
        print(text + "\n")
        f.write("\n")

def get_subtype_string(subtypes, API = False):
    '''returns a prefix for file headers or apa calls if API'''
    
    if len(subtypes)==1:
        return subtypes[0]
    
    if API:
        gap_char = "," 
    else:
        gap_char = "_"

    header = subtype[-1]
    for subtype in reversed(subtypes[0:-1]):
        header = f"{subtype}{gap_char}" + header
    return header 

def count_sequences(records):
    '''takes a Biopython FastaIO records and count sequences'''
    n = 0
    for _ in records:
        n+=1
    return n

def convert_string_to_seq(fasta_string):
    '''converts fasta text to a StringIO object'''
    fasta_io = StringIO(fasta_string) 
    records = SeqIO.parse(fasta_io, "fasta")
    return records

def clean(sequences, invalid_seq_path, cleaned_path):
    '''open output file and continously write valid sequences to output_fasta, 
    invalid sequences are written to invalid_seq_path for observation purposes'''
    ambiguous_aa = ["B", "J", "O", "U", "Z"]
    
    def is_valid_sequence(seq):
        '''Helper method. 
        return True if ambiguous_aa chars are in the string, 
        True if string is correctly formatted'''
        if seq[0]!="M": # sequence must start with M
            return False

        for char in ambiguous_aa:
            if char in seq:
                return False
        return True # if the loop completes, no invalid characters were found in inputted sequence

    invalid_sequence_counter=0
    with open (invalid_seq_path, "w+") as output_invalid: 
        with open(cleaned_path, "w+") as output:
            # Using the Biopython fasta parse we can read our fasta input
            for sequence in sequences:
                sequence_string = str(sequence.seq).upper()
                #write valid sequences to output
                if is_valid_sequence(sequence_string):
                    output.write(">" + str(sequence.description) + "\n" + str(sequence_string + "\n\n"))
                else:
                    invalid_sequence_counter+=1
                    output_invalid.write(">" + str(sequence.description) + "\n" + str(sequence_string + "\n\n"))
    log([f"{invalid_sequence_counter} invalid sequences were found"])


#### SHELL COMMANDS, CD-HIT, MUSCLE, PROTEIN-CONSENSUS-SEQUENCE #############
def run_shell(command):
    '''runs any given set of arguments on command line'''
    #print(args)
    with open(log_file, 'a') as f:
        process = subprocess.run(command, stdout=f, stderr = f, shell = True)

def run_shell(command):
    '''runs any given set of arguments on command line'''
    #print(args)
    with open(log_file, 'a') as f: # log file is defined at runtime
        process = subprocess.run(command, stdout=f, stderr = f, shell = True)
    
def cluster(cleaned_path, clustered_path, thresh = 0.99):
    '''performs cd-hit clustering'''
    log(["PERFORMING CD-HIT CLUSTERING", "-------------------------------------------------------"])
    terminal_command = f"{cd_hit_path} -i {cleaned_path} -o {clustered_path} -c {thresh}"
    run_shell(terminal_command)   

def align(clustered_path, aln_path):
    '''performs muscle alignment. May take a long time'''
    log(["MUSCLE ALIGNMENT", "-------------------------------------------------------"])
    terminal_command = f"muscle -align {clustered_path} -output {aln_path}"
    try:
        run_shell(terminal_command)
    except:
        print("ERROR: ensure muscle has been downloaded and activated as an executable using chmod")

def get_consensus(aln_path, consensus_path):
    '''CONSENSUS SEQUENCE AND ENTROPY ANALYSIS - MODIFIED PROTEIN-CONSENSUS-SEQUENCE PROGRAM

    runs protein-consensus-sequence program from https://github.com/msternke/protein-consensus-sequence.git
   
    NOTE: I have included a modified version of this program in the vaccine_design package,
            no need to download any additional software (besides cd-hit and muscle)
   
    
    Outputs
    ---------------------------------------
    1) calculates consensus sequence with gap removal
    2) generate residues frequency
    3) Entropy data at each position
    4) Entropy plots 
    '''
    terminal_command = f"{consensus_path}/consensus.py -i {aln_path} -o consensus.fasta -c 1"
    try:
        run_shell(terminal_command)
    except:
        log(["ensure that you have downloaded protein-consensus-sequence and it has been placedd in working directory"])

def process_consensus(p, s, d):
    '''get residue frequencies specifially for consensus amino acid and save to a file'''
    seq_record = SeqIO.parse(s, "fasta")
    first_sequence = (next(seq_record))
    sequence =  first_sequence.seq.split("MSA")[0]
    #print(sequence)
    description =  first_sequence.description

    # write sequence to file
    with open (d, "w") as f:
        f.write(">" + str(description) + "\n")
        f.write(str(sequence))

    df = pd.read_csv(p)
    frequency_dict={"pos":[],
                "aa":[],
               "freq":[]}
    for idx,char in enumerate(list(sequence)):
        freq = df.loc[idx,char]
        pos = idx+1
        frequency_dict["pos"].append(pos)
        frequency_dict["aa"].append(char)
        frequency_dict["freq"].append(freq)

    # convert data to dataframe and save as csv
    output_df = pd.DataFrame(frequency_dict)
    output_df.set_index("pos",inplace=True)
    save_freq_chart_head = os.path.split(p)[0]
    save_freq_chart_path = os.path.join(save_freq_chart_head, "consensus_info.csv")
    output_df.to_csv(save_freq_chart_path)



def main(subtypes):
    # Step 0
    # Set up
    parent_path = os.getcwd() # returns absolute path to this directory

    subtypes.sort()
    subtype_string = get_subtype_string(subtypes) # generate an file header to identify this session

    # define file paths for outputs
    consensus_path =  os.path.join(parent_path, "protein-consensus-sequence")
    project_dir = os.path.join(parent_path, f"data/{subtype_string}") # location of session dir
    data_path = os.path.join(parent_path, "data")
    clustered_path = "clustered.fasta" 
    invalid_seq_path = "invalid.fasta"
    cleaned_path = "seqs.fasta"
    aln_path = "aln.fasta"
    consensus_out_path = os.path.join(project_dir, "Outputs")

    # create project directories and log file about sequence data
    os.makedirs(project_dir, exist_ok=True)
    
    # change current working directory
    os.chdir(project_dir)

    # begin logging to log.txt
    log([subtype_string], newline = False)

    # Step 1 
    # Get data from IRD
    # API Call using requests
    subtype_string_api = get_subtype_string(subtypes, API = True)

    base_url = 'https://www.fludb.org/brc/api/sequence'

    # must divide search into pre 2008 and post 2008 to avoid overload
    # 10,000 sequences max can be quiered at a time
    url1 = f"?datatype=protein&completeseq=y&host=human&family=influenza&toyear=2008&flutype=A&protein=HA&subtype={subtype_string_api}&metadata=uniprotAcc,strainName,subtype&output=fasta"
    url2 = f"?datatype=protein&completeseq=y&host=human&family=influenza&fromyear=2009&flutype=A&protein=HA&subtype={subtype_string_api}&metadata=uniprotAcc,strainName,subtype&output=fasta"

    # make API call
    response1 = requests.get(base_url+url1)
    response2 = requests.get(base_url+url2)
    log(["IRDB API QUERY", response1.url, response2.url])

    sequence_string = response1.text + "\n" + response2.text

    # first, convert sequences to Biopython records object
    sequences = convert_string_to_seq(sequence_string)

    # count sequences and write to file
    n = count_sequences(sequences)
    log([f"{n} sequences were quiered from the IRD"])

    # perform cleaning using function above
    # After iterating through Biopython FASTAIO generator, it's empty
    # must restore the sequence record
    sequences = convert_string_to_seq(sequence_string)
    clean(sequences, invalid_seq_path, cleaned_path)

    # cluster sequences
    cluster(cleaned_path, clustered_path)

    # count clusters for log file 
    clustered_records = SeqIO.parse(clustered_path, "fasta")
    num_clusters = count_sequences(clustered_records)
    log([f"{num_clusters} clusters were created"])

    # align sequences
    # output of muscle is sent to log file
    # BUG: for some reason, output is not written to file
    # BUG SOLVED: muscle uses stderr to print output, not stdout like other command line tools
    align(clustered_path, aln_path)

    # get the consensus sequence
    # entropy is calculated 
    get_consensus(aln_path, consensus_path)

    # get residue frequencies
    freq_chart_path = os.path.join(consensus_out_path, "consensus_residueFrequencies.csv")
    consensus_seq_path =  os.path.join(consensus_out_path, "consensus_consensus_output.txt")
    
    # move entropy data to data directory
    entropy_data_path = os.path.join(data_path, f"{subtype_string}_entropy.csv")
    new_ent_path = shutil.copy("entropy.csv", data_path)
    os.rename(new_ent_path, entropy_data_path)
    
    # write data to data directory
    consensus_seq_data_path = os.path.join(data_path, f"{subtype_string}_consensus.fasta")
    process_consensus(freq_chart_path, consensus_seq_path, consensus_seq_data_path)

    # in case the user wants to loop over subtypes, set working directory back to vaccine_design path
    os.chdir(parent_path)

if __name__=="__main__":
    # must pass subtypes as list
    # change this for your job
    subtypes_to_run = [
                        ["H7"]
    ]

    for subtype in subtypes_to_run:
        print("running consensus sequence generation for", subtype)
        main(subtype)

