from io import StringIO
from Bio import SeqIO
from pymol import cmd
import sys
import os
from utils import get_start_and_end
from pymol_utils import remove_non_RBS


def align_all_RBS_to_all(object_list=None,selection='name ca',cutoff=2,
                        cycles=5,debug=0,full_matrix=1,method='align', 
                        RBS_path = "data/consensus_RBS.fasta" ):
  """
  Berk Yalcinkaya
  Adapted from Robert L. Campbell (rlc1@queensu.ca)
  Feel free to do whatever you like with this code.

  Isolates the RBS of each HA protein using pairwise alignment
  and Structurally Aligns all RBS's to all others in the given selection

  usage:
    align_all_RBS_to_all [object_list][selection][cutoff=2][cycles=5][debug=0][full_matrix=0][method='align']

        where method can be align, super or cealign

        where selection can be used to specify a particular selection
        (e.g. one chain of each protein),

        cutoff and cycles are options passed to the align or super command.
    
    First, RBS is isolated based on given RBS sequence (RBS path) via pairwise alingment.

    By default, all objects are aligned to all others using just the C-alpha atoms.
    You can specify a list of objects to use instead with the object_list option.

    Setting debug=1 prints more information to the terminal or external GUI.

    Setting full_matrix=1 prints out the full symmetric matrix, rather than
    simply the top-half matrix


    Example:
      align_all_to_all object_list=name1 name2 name3 name4, selection=c. a & n. ca, full_matrix=1
  """

  # MUST BE LOCATED IN VACCINE_DESIGN DIR TO USE
  RBS_seq_file = RBS_path
  RBS_record = SeqIO.parse(RBS_seq_file, "fasta")
  RBS_seq = (next(RBS_record)).seq    

  cutoff = int(cutoff)
  full_matrix = int(full_matrix)
  cycles = int(cycles)
  debug=int(debug)

  if not object_list:
    object_list = cmd.get_names()
  else:
    object_list = object_list.replace('[','').replace(']','').replace(',',' ').split()
  
  #isolate RBS for each object to align
  for object in object_list:
    pymol_fasta_string = cmd.get_fastastr(object)
    fasta_io = StringIO(pymol_fasta_string) 
    sample_record = next(SeqIO.parse(fasta_io, "fasta"))
    sample_seq = sample_record.seq
    start,end = get_start_and_end(sample_seq, RBS_seq)
    remove_non_RBS(object,start,end)

  rmsd = {}
  rmsd_list = []
  #  print object_list
  for i in range(len(object_list)):
    for j in range(i+1,len(object_list)):
      if method == 'align':
        rms = cmd.align('%s & %s' % (object_list[j],selection),'%s & %s' % (object_list[i],selection),cutoff=cutoff,cycles=cycles)
      elif method == 'super':
        rms = cmd.super('%s & %s' % (object_list[j],selection),'%s & %s' % (object_list[i],selection),cutoff=cutoff,cycles=cycles)
      elif method == 'cealign':
        rmsdict = cmd.cealign('%s & %s' % (object_list[i],selection),'%s & %s' % (object_list[j],selection))
        rms = [rmsdict['RMSD'],rmsdict['alignment_length'],1,0,0]
      elif method == 'rms_cur':
#        print'mobile: %s & %s' % (object_list[j],selection),'target: %s & %s' % (object_list[i],selection)
        num_atoms = cmd.select('junkselection','%s & %s' % (object_list[j],selection))
        cmd.delete('junkselection')
        rms = [cmd.rms_cur('%s & %s' % (object_list[j],selection),'%s & %s' % (object_list[i],selection)),num_atoms]

      else:
        print("Method: ",method)
        print("only 'align', 'super', 'cealign' and 'rms_cur' are accepted as methods")
        sys.exit(-1)

      rmsd.setdefault(object_list[i],{})[object_list[j]] = rms[0]
      rmsd_list.append((object_list[j],object_list[i],rms[0],rms[1]))
      if debug and method != 'rms_cur':
        print("Alignment of %s to %s:" % (object_list[j],object_list[i]))
        print("     Initial RMS: %6.3f for %d atoms" % (rms[3],rms[4]))
        print("     Final RMS: %6.3f for %d atoms after %d cycles\n" % (rms[0],rms[1],rms[2]))

#  sort RMSD list
  rmsd_list.sort(key=lambda x: x[2])

# loop over dictionary and print out matrix of final rms values
  if debug:
    for object_name in object_list[:-1]:
      print("%s: %s" % (object_name,str(rmsd[object_name])))

    print("\nSorted from best match to worst:")
    for r in rmsd_list:
      print("%s to %s: %6.3f using %d atoms" % r)
    print("")

  print("%6s" % " ", end=' ')
  if full_matrix:
# fill in other half of matrix
    for i in range(len(object_list)):
      for j in range(i+1,len(object_list)):
        rmsd.setdefault(object_list[j],{})[object_list[i]] = rmsd[object_list[i]][object_list[j]]
      rmsd[object_list[i]][object_list[i]] = 0

    for i in range(len(rmsd)):
      print("%6s" % object_list[i], end=' ')
    print("")
    for i in range(len(object_list)):
      print("%6s" % object_list[i], end=' ')
      for j in range(len(object_list)):
        print("%6.3f" % (rmsd[object_list[i]][object_list[j]]), end=' ')
      print("")
  else:
    
    with open('test.txt','w') as file_out:
      for i in range(len(rmsd)):
        print("%6s" % object_list[i+1], end=' ')
        file_out.write(f"{object_list[i+1]:6s}")
      print("")
      file_out.write(f"\n")

      for i in range(len(object_list)):
        print("%6s" % object_list[i], end=' ')
        file_out.write(f"{object_list[i]:6s}")
        for k in range(i):
          print("%6s" % " ", end=' ')
          file_out.write(f"{' ':6s}")
        for j in range(i+1,len(object_list)):
          print("%6.3f" % (rmsd[object_list[i]][object_list[j]]), end=' ')
          file_out.write(f"{rmsd[object_list[i]][object_list[j]]:6.3f}")
        print("")
        file_out.write(f"\n")

    print(f"outfile:{os.path.join(os.getcwd(),'test.txt')}")

cmd.extend('align_all_RBS_to_all',align_all_RBS_to_all)
        