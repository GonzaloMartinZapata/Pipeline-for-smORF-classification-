import os 
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
import argparse
parser= argparse.ArgumentParser ()
parser.add_argument ("--folder", "-f" ,dest= "folder", required= True, help= "folder with the genomes" , nargs = 1)                 
args= parser.parse_args ()

def annotated_smorfs (genomes, out_files):
    with open (out_files, "a") as out_handle:
        print (f"processing {genomes}")
        for record in SeqIO.parse (genomes , "gb") :
            for feature in record.features :
                if feature.type == "CDS" :
                    if "locus_tag" in feature.qualifiers.keys () :
                        locus_tag = feature.qualifiers ["locus_tag"] [0]

                    if "translation" in feature.qualifiers.keys () :
                        translation = feature.qualifiers ["translation"] [0]

                    if "translation" in feature.qualifiers.keys () and len (feature.qualifiers ["translation"] [0]) <= 70 :
                        out_handle.write (">" + locus_tag + "\n" + translation + "\n")
                        

path= os.getcwd ()
folder= args.folder [0]
working_directory= os.path.join (path, folder)
results= "Annotated_smORFs"
results_directory= os.path.join (working_directory, results)
if not os.path.exists (results_directory):
    os.mkdir (results_directory)

files = glob.glob (os.path.join (working_directory, "*gbff"))
for genomes in files:
    basename, extension= os.path.splitext (os.path.basename (genomes))
    txt_files= basename + ".fa"
    out_files= os.path.join (results_directory, txt_files)
    annotated_smorfs (genomes, out_files)
    