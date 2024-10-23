import pandas as pd
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings ("ignore")
import argparse
parser= argparse.ArgumentParser ()
parser.add_argument ("--results", "-r", dest= "results", required=True, help= "path to the file with all the smORFs of a Bachteria", nargs=1)
parser.add_argument ("--genomes", "-g", dest="genomes", required=True, type=int, help="number of genomes where the smORFs were extracted", nargs=1 )
args= parser.parse_args ()

results= args.results [0]

argument, extension= os.path.splitext (results)

df=pd.read_csv (results)

priority_order = ["Conserved", "Core", "Soft core", "Shell", "Cloud", "Probably not conserved"]

#Average of smORFs in a genome
smORFS_x_genome= (len (df)/args.genomes [0])
smORFS_x_genome= round (smORFS_x_genome)

df["Results_column"] = pd.Categorical(df["Results_column"], categories=priority_order, ordered=True)

genes50 = df[df["Results_column"] < "Shell"]

conserved = df[df["Results_column"] == "Conserved"]

core= df[df["Results_column"] == "Core"]

soft_core= df[df["Results_column"] == "Soft core"]

#Calculate the percentage of smORFs present in 50% of the genomes
percentage_smorfs= (len (genes50)*100)/(len (df) ) 
percentage_smorfs= round (percentage_smorfs, 2)

#Calculate the average smORFs50 in a genome
smorfs50_x_genome= (len (genes50)/ args.genomes [0]) 
smorfs50_x_genome= round (smorfs50_x_genome)

#Calculate the average of conserved,core and soft core genes in a single genome
conserved_x_genome= (len (conserved)/ args.genomes [0])
conserved_x_genome= round (conserved_x_genome)

core_x_genome= (len (core)/ args.genomes [0])
core_x_genome= round (core_x_genome)

soft_core_x_genome= (len (soft_core)/ args.genomes [0])
soft_core_x_genome= round (soft_core_x_genome)

#Number of probably conserved genes that are annotated as "hypothetical Protein"
hypothetical_50= genes50 [genes50["product"]== "hypothetical protein"]

hypothetical_conserved= conserved [conserved["product"]== "hypothetical protein"]

hypothetical_core= core [core["product"]== "hypothetical protein"]

hypothetical_soft_core= soft_core [soft_core["product"]== "hypothetical protein"]

#smORFs classified in the same category in the three programs
equal_class = df[(df["Classification_of"] == df["Classification_panaroo"]) & (df["Classification_of"] == df["Classification_ognc"])]

percentage_equal_class= (len (equal_class)*100)/(len (df) ) 
percentage_equal_class= round (percentage_equal_class, 2)


#The name of the file has to be the same as the name of the organism
argument, extension= os.path.splitext (args.results [0])

results_file= argument + "_results" + ".txt"

#Writing of the results file
with open (results_file, "a") as file:
    file.write ("In " + str(argument) + "," + str (len (df)) + " smORFs have been found. \n")
    file.write ("Out of them, " + str (len (genes50)) + " are present in at least 50% of the genomes. (" + str(percentage_smorfs) + "%) \n")
    file.write ("In average, there are " + str (smORFS_x_genome) + " smORFs in a genome, out of which " + str(smorfs50_x_genome) + " are probably conserved. \n")
    file.write (str(conserved_x_genome)+ " are conserved, " + str(core_x_genome) + " are core, and " + str (soft_core_x_genome) + " are soft core. \n")
    file.write("Out of the " + str(len(genes50)) + " smORFs present in at least 50% of the genomes, " + str(len(hypothetical_50)) + " are annotated as hypothetical proteins. \n")
    file.write ("(" + str(len (hypothetical_conserved)) + " conserved, " + str (len (hypothetical_core)) + " core, and " + str (len (hypothetical_soft_core)) + " soft core). \n")
    file.write (str (len (equal_class)) + " smORFs, were classified in the same group by the three programs. (" + str(percentage_equal_class) + "%) \n")

def piechart (df):

    # Crea una nueva columna con el orden específico
    df['Results_column'] = pd.Categorical(df['Results_column'], categories=priority_order, ordered=True)

    classification= df.groupby ("Results_column") ["genome"].count ().reset_index ()
    frecuencies= classification.iloc [:,1]
    category= classification.iloc [:, 0]
    colors = ["deepskyblue", "lightblue", "lightgreen", "violet", "lightcoral", "crimson"]

    pie, _ = plt.pie(frecuencies, startangle=90, colors=colors)

    porcent= (classification["genome"]/classification ["genome"].sum ()*100)
    labels = ['{0} - {1:1.2f} %'.format(i,j) for i,j in zip(category, porcent)]

    plt.legend(pie, labels, loc="upper left", bbox_to_anchor= (1,1))
    plt.savefig (argument + "_smORFs_classified" + ".svg", bbox_inches= "tight", format= "svg")


piechart (df)