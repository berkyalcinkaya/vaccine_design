from pymol import cmd
import pandas as pd

'''in PyMOL, use by:
run color_by_identity.py'''

print("running script")

def color_by_conservation():

    def get_color_name(freq):
        freq = freq*100
        if freq<=50:
            return "color_0"
        if freq <= 60:
            return "color_1"
        if freq<= 70:
            return "color_2"
        if freq<= 80:
            return "color_3"
        if freq<= 85:
            return "color_4"
        if freq<=90:
            return "color_5"
        if freq<=95:
            return "color_6"
        if freq<=98:
            return "color_7"
        if freq<=100:
            return "color_8"

    def get_sample_res(ref_idx):
        idx=-1
        for char in sample_seq[0:ref_idx+1]:
            if char!="-":
                idx+=1
        if idx == -1:
            raise KeyError
        else:
            return idx


    # rainbow colors should be:
        # gray, light blue, blue, teal, green, yellow, orange, red, magenta
        # to correspond to identity of <= :
        # 50, 60, 70, 80, 85, 90, 95, 98, 100
    #freq_vals = [50, 60, 70, 80, 85, 90, 95, 98, 100]
    colors={}
    colors["rainbow"] = [[0.5, 0.5, 0.5] , [0.38,0.38,0.72] , [0,0.30,0.70], [0.06,0.68,0.40], 
                [0.21,0.66,0], [0.85,0.66,0], [0.94,0.32,0], [0.88,0,0], [1,0,0.55]]
    colors["green"] = [ [0.63,0.63,0.63] , [0.50,0.68,0.56] , [0.42,0.71,0.53], 
				   [0.35,0.74,0.49] , [0.26,0.77,0.44] , [0.19,0.80,0.41], 
				   [0.12,0.83,0.37] , [0.01,0.87,0.31] , [0,1,0.83] ]

    for i,color in enumerate(colors["green"]):
        color_string = f"color_{i}"
        cmd.set_color(color_string,color)


    sample_seq = "---------------GSDTICIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDSHNGKLCRLKGIAPLQLGKCNIAGWILGNPECESLLSKRSWSYIAETPNSENGTCYPGDFADYEELREQLSSVSSFERFEIFPKERSWPKHNITRGVTAACSHAGKSSFYKNLLWLTETNGSYPKLSKSYVNNKEKEVLVLWGVHHPSNIEDQKTLYRKENAYVSVVSSNYNRRFTPEIAERPKVRGQAGRMNYYWTLLEPGDTIIFEANGNLIAPWYAFALSRGFGSGIITSNASMDECDTKCQTPQGAINSSLPFQNIHPVTIGECPKYVKSTKLRMVTGLRNIPSIQSRGLFGAIAGFIEGGWTGMIDGWYGYHHQNEQGSGYAADQKSTQNAINGITNKVNSVIEKMNTQFTAVGKEFNKLEKRMENLNKKVDDGFLDIWTYNAELLVLLENERTLDFHDSNVKNLYEKVKNQLRNNAKEIGNGCFEFYHKCNNECMESVKNGTYDYPKYSEEFLVPR----------------------------------------------------"
    sample_seq =  sample_seq.replace("\n",  "").replace(" ", "").strip()


    ref_seq = '''MKTILAVLLYIFCLVSADTLCIGYHANNSTDTVDTILEKNVTVTHSVNLLEDSHNGKLCD
    LPGVAPLDLGNCTLAGWLLGNPECDGLQNAKSWSYIVERSKADNGTCYPGDFPDYEELRE
    QLSSVSSFERFEIFPKSSSWPNHT-QNGVSAACPHAGAKSFYSNLNWLTKKGNSYPALNV
    TYPNNKGKEVLVLWGVHHPSTDADQQSLYQNASGYVTVSTSRYSQKFIPEIASRPKVRDQ
    EGRINYYWTLVKPGDIILFEATGNLIAPRYAFKIERSGGSGIIRSDAPIGKCNTECQTPN
    GAINNSLPFQNVHPITIGACPKYVKSTKLRLATGLRNVPSIQSRGLFGAIAGFIEGGWTG
    MVDGWYGYHHQNEQGSGYAADLKSTQNAIDQITNKVNSVIEKMNTQFTAVGKEFNELEGR
    IENLNKKVDDGFLDIWTYNAELLVLLENERTLDLHDSNVKNLYEKVRKQLRNNAKEIGNG
    CFEFYHKCDNACMESVRNGTYDYPKYSEEAKLNREQIDGVKLESGYKYQILAIYSTVASS
    LVLCVSLGAISFWMCSNGSLQCRICI'''
    ref_seq = ref_seq.replace("\n",  "").replace(" ", "").strip()

    # load frequnecy table
    consensus_data = pd.read_csv("/Users/berk/vaccine/data/other_ref_sequences/muscle_align/H1_H3_H5/consensus_info.csv", index_col="pos")

    ref_aa_pos = 0 # because of gaps in ref seq, idx is not the same as the aa pos
    for idx,aa in enumerate(ref_seq):
        
        # skip over gaps in reference sequence
        if aa=="-":
            continue
        
        # if the two sequences have the same aa, perform coloring based on aa frequency at that position
        elif sample_seq[idx]==aa:
            sample_idx = get_sample_res(idx)
            sample_aa_pos = sample_idx+1 # aa pos is 1 indexed
            
            freq = consensus_data.iloc[ref_aa_pos]["freq"]
                
            #print(i+1, sample_seq[i],aa, round(freq,3))
            #print(i+1,samepl_aa_pos,freq)
            
            color_name = get_color_name(freq)
            selection_string = f"resi {sample_aa_pos}" 
            cmd.select("selection",selection_string)
            cmd.color(color_name, "selection")
            
            # increment reference sequence aa position only if aa is not a gap char
            ref_aa_pos+=1

cmd.extend("color_by_conservation", color_by_conservation)
color_by_conservation()