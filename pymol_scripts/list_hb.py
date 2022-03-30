from itertools import combinations
import pandas as pd
import os

from pymol import cmd
from pymol import stored


print("running script")

def list_hb(selection, selection2=None, cutoff=3.2,
        angle=45, mode=1, hb_list_name='hbonds',print_distances=1,write_distances_file=None):
  """
  USAGE

  list_hb selection, [selection2 (default=None)], [cutoff (default=3.2)],
                     [angle (default=45)], [mode (default=1)],
                     [hb_list_name (default='hbonds')], [print_distances (default=1)]
                     [write_distances (default=None)]

  The script automatically adds a requirement that atoms in the
  selection (and selection2 if used) must be either of the elements N or
  O.

  If mode is set to 0 instead of the default value 1, then no angle
  cutoff is used, otherwise the angle cutoff is used and defaults to 45
  degrees.

  e.g.
  To get a list of all H-bonds within chain A of an object and to save it to a file
    list_hb 1abc & c. a &! r. hoh, cutoff=3.2, hb_list_name=abc-hbonds,write_distances_file=abc-hbonds.dat

  To get a list of H-bonds between chain B and everything else:
    list_hb 1tl9 & c. b, 1tl9 &! c. b
  """

  if write_distances_file:
    hb_data = open(write_distances_file, 'w')

  cutoff = float(cutoff)
  angle = float(angle)
  mode = int(mode)
  print_distances=int(print_distances)
  selection = str(selection)
  selection2 = str(selection2)

  #define hydrogen bond donors and accpetoprs
  acceptors = ['(e. C & ! sc. )', '(r. asn & n. OD1)', '(r. asp & n. OD1+OD2)',
  '(r. glu & n. OE1+OE2)', '(r. his & n. ND1)', '(r. met & n. SD)', '(r. ser & n. OG)',
  '(r. thr & n. OG1)', '(r. tyr & n. OH)']
  accept_str = " or ".join(acceptors)
  acceptor = " & ({})".format(accept_str)

  donors = ['(n. N & ! sc. & nbr. hydro)', '(r. arg & n. NE+NH1+NH2)', '(r. asn & n. ND2)', '(r. cys & n. SG)',
  '(r. gln & n. NE2)', '(r. his & n. NE2+ND1)', '(r. ser & n. OG)', '(r. lys & n. NZ)', '(r. thr &n. OG1)',
  '(r. trp & n. NE1)', '(r. tyr & n. OH)']
  don_str = " or ".join(donors)
  donor = " & ({})".format(don_str)


  #get pairs
  if selection2==None:
      hb = cmd.find_pairs(selection+donor, selection+acceptor, mode=mode,
              cutoff=cutoff, angle=angle)
  else:
      hb_first = cmd.find_pairs(selection+donor, selection2+acceptor, mode=mode,
              cutoff=cutoff, angle=angle)
      hb_second = cmd.find_pairs(selection+acceptor, selection2+donor, mode=mode,
              cutoff=cutoff, angle=angle)
      hb = hb_first + hb_second

  #report empty list
  if len(hb) == 0:
    print("No pairs detected")
    print("_________________")


  # convert hb list to set to remove duplicates
  hb_set = set()
  for atoms in hb:
    a = [atoms[0],atoms[1]]
    a.sort()
    hb_set.add(tuple(a))

# convert set back to list and sort for easier reading
  hb = list(hb_set)
  hb.sort(key=lambda x: x[0][1])

  stored.listA = []
  stored.listB = []
  stored.listC = []

  #get and store attributes of atoms in interaction pairs
  for pairs in hb:
    cmd.iterate("%s and index %s" % (pairs[0][0], pairs[0][1]),
            'stored.listA.append( "%1s/%3s`%s/%s/%i " % (chain, resn, resi, name, ID),)')

    cmd.iterate("%s and index %s" % (pairs[1][0], pairs[1][1]),
            'stored.listB.append( "%1s/%3s`%s/%s/%i " % (chain, resn, resi, name, ID),)')

    stored.listC.append(cmd.distance(hb_list_name, "%s and index %s" %
                (pairs[0][0], pairs[0][1]), "%s and index %s" % (pairs[1][0],
                    pairs[1][1])))

  #write to file,print, get interaction map
  for line in enumerate(stored.listA):
    if print_distances:
      print("%s   %s   %.2f" % (stored.listA[line[0]], stored.listB[line[0]], stored.listC[line[0]]))
    if write_distances_file:
      hb_data.write("%s   %s   %.2f\n" % (stored.listA[line[0]], stored.listB[line[0]], stored.listC[line[0]]))
  if write_distances_file:
    hb_data.close()
    generateMap(write_distances_file) #write data to csv file as map

  return stored.listA, stored.listB, stored.listC


#list_hb_btwn chains
def list_hb_btwn_chains(selection, cutoff=3.2,
        angle=45, mode=1, hb_list_name='hbonds',print_distances=1,write_distances_file=None):
  """
  USAGE

  see list_hb
  input:
  selection - protein whoses interchain interactions are to be analyzed
  cutoff - max distance between atoms
  angle - allowed difference from 180
  """
  if write_distances_file:
    hb_data = open(write_distances_file, 'w')

  cutoff = float(cutoff)
  angle = float(angle)
  mode = int(mode)
  print_distances=int(print_distances)
  selection = str(selection)

  #define hydrogen bond donors and accpetoprs
  acceptors = ['(e. O)','(e. C & ! sc. )', '(r. asn & n. OD1)', '(r. asp & n. OD1+OD2)', '(r. glu & n. OE1+OE2)',
  '(r. his & n. ND1)', '(r. met & n. SD)', '(r. ser & n. OG)', '(r. thr & n. OG1)', '(r. tyr & n. OH)']
  accept_str = " or ".join(acceptors)
  acceptor = " & ({})".format(accept_str)

  donors = ['(n. N & ! sc. & nbr. hydro)', '(r. arg & n. NE+NH1+NH2)', '(r. asn & n. ND2)', '(r. cys & n. SG)',
  '(r. gln & n. NE2)', '(r. his & n. NE2+ND1)', '(r. ser & n. OG)', '(r. lys & n. NZ)',
  '(r. thr & n. OG1)', '(r. trp & n. NE1)', '(r. tyr & n. OH)' ]
  don_str = " or ".join(donors)
  donor = " & ({})".format(don_str)

  #generate chain pairs to be analyzed
  chains_lst = cmd.get_chains(selection)
  chain_pairs = list(combinations(chains_lst, 2))
  print(chain_pairs)

  #find interaction atom pairs from chain pairs
  hb=[]
  for pair in chain_pairs:

      selection1 = selection + " & c. " + str(pair[0])
      selection2 = selection + " & c. " + str(pair[1])

      hb_first = cmd.find_pairs(selection1+donor, selection2+acceptor, mode=mode,
              cutoff=cutoff, angle=angle)
      hb_second = cmd.find_pairs(selection2+donor, selection1+acceptor, mode=mode,
              cutoff=cutoff, angle=angle)
      hb = hb + hb_first + hb_second

  # convert hb list to set to remove duplicates
  hb_set = set()
  for atoms in hb:
    a = [atoms[0],atoms[1]]
    a.sort()
    hb_set.add(tuple(a))

  #convert set back to list and sort for easier reading
  hb = list(hb_set)
  hb.sort(key=lambda x: x[0][1])

  stored.listA = []
  stored.listB = []
  stored.listC = []

  #get atom pair attributes from PyMOL
  for pairs in hb:
    cmd.iterate("%s and index %s" % (pairs[0][0], pairs[0][1]),
            'stored.listA.append( "%1s/%3s`%s/%s/%i " % (chain, resn, resi, name, ID),)')

    cmd.iterate("%s and index %s" % (pairs[1][0], pairs[1][1]),
            'stored.listB.append( "%1s/%3s`%s/%s/%i " % (chain, resn, resi, name, ID),)')

    stored.listC.append(cmd.distance(hb_list_name, "%s and index %s" %
                (pairs[0][0], pairs[0][1]), "%s and index %s" % (pairs[1][0],
                    pairs[1][1])))

  #process interaction data
  for line in enumerate(stored.listA):
    if print_distances:
      print("%s   %s   %.2f" % (stored.listA[line[0]], stored.listB[line[0]], stored.listC[line[0]]))
    if write_distances_file:
      hb_data.write("%s   %s   %.2f\n" % (stored.listA[line[0]], stored.listB[line[0]], stored.listC[line[0]]))
  if write_distances_file:
    hb_data.close()
    generateMap(write_distances_file) # generate interaction map

  return stored.listA, stored.listB, stored.listC


def generateMap(path):
  '''
  opens the file created above and parses the text to create a cleaner interaction map

  input: path to output text file
  output: none
  '''
  output_filepath = os.path.splitext(path)[0]+'.csv'

  with open(path, 'r') as text:
    lines_lst = text.readlines()

  lst = []
  for line in lines_lst[1:]:
    #print(line)
    components = line.split("    ")
    chain_1 = components[0].split("/")[0]
    chain_1_pos = components[0].split("/")[1].split("`")[1]
    chain_1_aa = components[0].split("/")[1].split("`")[0]
    chain_2 = components[1].split("/")[0]
    chain_2_pos = components[1].split("/")[1].split("`")[1]
    chain_2_aa = components[1].split("/")[1].split("`")[0]
    #dist = components[2]
    dict = {"Chain": chain_1,
           "Chain Position":int(chain_1_pos),
           "Chain Res":chain_1_aa,
           "Interacting Chain": chain_2,
           "Interacting Chain Position":int(chain_2_pos),
           "Interacting Chain Res":chain_2_aa}
    if dict not in lst:
        #print(dict)
        lst.append(dict)
  df = pd.DataFrame(lst)
  df = df.sort_values(by = ["Chain", "Chain Position"], ascending=True)
  df.to_csv(output_filepath, index = False)

  return None


#run extend commands
cmd.extend("list_hb", list_hb)
print("The function list_hb has been added to your PyMOL PATH.")
cmd.extend("list_hb_btwn_chains", list_hb_btwn_chains)
print("The function list_hb_btwn_chains has been added to your PyMOL PATH.")
