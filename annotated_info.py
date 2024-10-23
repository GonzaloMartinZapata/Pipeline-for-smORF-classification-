import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--folder", "-f", dest="folder", required=True, help="folder with the genomes")
args = parser.parse_args()

def annotated_info(genomes):
    print(f"processing {genomes}")
    tables = pd.DataFrame(columns=["genome", "locus_tag", "prot_seq", "len", "product", "function", "note"])
    basename, extension = os.path.splitext(os.path.basename(genomes))
    for record in SeqIO.parse(genomes, "gb"):
        for feature in record.features:
            if "locus_tag" in feature.qualifiers.keys():
                locus_tag = feature.qualifiers["locus_tag"][0]
                    
            if "translation" in feature.qualifiers.keys():
                translation = feature.qualifiers["translation"][0]
            if not "translation" in feature.qualifiers.keys():
                translation = " "

            if "function" in feature.qualifiers.keys():
                function = feature.qualifiers["function"][0]
            if not "function" in feature.qualifiers.keys():
                function = " "    
                    
            if "product" in feature.qualifiers.keys():
                product = feature.qualifiers["product"][0]
            if not "product" in feature.qualifiers.keys():
                product = " "        
            if "note" in feature.qualifiers.keys():
                note = feature.qualifiers["note"][0]
            if not "note" in feature.qualifiers.keys(): 
                note = " "

            if "translation" in feature.qualifiers.keys() and (len(feature.qualifiers["translation"][0])) <= 70:
                dict = {"genome": [basename], "locus_tag": [locus_tag], "prot_seq": [translation], "len": [str(len(translation))], "product": [product], "function": [function], "note": [note]}
                df = pd.DataFrame(dict)

                tables = pd.concat([tables, df])

    tables = tables.reset_index(drop=True)
    return tables 

def functions(data):
    if data["product"].empty:
        print("No hay datos en la columna 'product' para procesar.")
        return
    
    top_categories = data["product"].value_counts().head(10)
    other_categories = data["product"].value_counts().tail(-10).sum()
    
    if top_categories.empty:
        print("No hay categorías para graficar.")
        return
    
    df_combined_categories = pd.concat([top_categories, pd.Series({"Others": other_categories})])
    
    if df_combined_categories.empty:
        print("No hay datos combinados para graficar.")
        return

    colores = ['lightcoral', 'lightskyblue', 'lightgreen', 'lightpink', 'lightsalmon', 'lightyellow', 'lightblue', 'lightcyan', 'lightgray', 'orange', 'lightseagreen']
    pie, _ = plt.pie(df_combined_categories, startangle=90, colors=colores)

    porcentajes = (df_combined_categories / df_combined_categories.sum() * 100)
    labels = [f'{categoria} - {porcentaje:.1f}%' for categoria, porcentaje in zip(df_combined_categories.index, porcentajes)]

    plt.legend(pie, labels, loc="upper left", bbox_to_anchor=(1, 1))

    plt.savefig("products_smorfs_" + args.folder + ".svg", bbox_inches="tight", format= "svg")
    plt.close()

def histogram(data):
    df = data
    len_column = df["len"]
    in_order = len_column.sort_values()
    range_bins = [10, 20, 30, 40, 50, 60, 70]

    sns.displot(len_column, bins=range_bins, kde=False)
    plt.xticks()
    plt.ylabel("Frecuencies")
    plt.xlabel("SEPs size")
    plt.savefig("histogram_SEPs_" + args.folder + ".svg", format= "svg")
    plt.close()

def boxplot(data):
    unique_genomes = data['genome'].unique()

    bins_list = []
    count_list = []

    for genome in unique_genomes:
        subset_df = data[data['genome'] == genome]

        subset_df["bins"] = pd.cut(subset_df["len"], bins=[0, 10, 20, 30, 40, 50, 60, 70])

        bins = subset_df.groupby("bins")["len"].count().reset_index()
        bins_list.extend(bins['bins'])
        count_list.extend(bins['len'])

    new_df = pd.DataFrame({'bins': bins_list, 'count': count_list})

    sns.boxplot(data=new_df, x="bins", y="count")
    plt.xlabel("SEPs size (aa)")
    plt.ylabel("Count")
    plt.xticks(rotation=45)
    plt.savefig("boxplot_smorfs_" + args.folder + ".svg", bbox_inches="tight", format= "svg")
    plt.close()

path = os.getcwd()
folder = args.folder
working_directory = os.path.join(path, folder)
results = "Annotated_info"
results_directory = os.path.join(working_directory, results)
if not os.path.exists(results_directory):
    os.mkdir(results_directory)

files = glob.glob(os.path.join(working_directory, "*gbff"))
combined_df = pd.DataFrame(columns=["genome", "locus_tag", "prot_seq", "len", "product", "function", "note"])
for genomes in files:
    tables = annotated_info(genomes)
    combined_df = pd.concat([combined_df, tables], ignore_index=True)
    combined_df["product"] = combined_df["product"].str.lower()

combined_df.to_csv(os.path.join(results_directory, args.folder + ".csv"), index=False)

os.chdir(results_directory)
csv_file = glob.glob(os.path.join(results_directory, "*.csv"))
for annotated_data in csv_file:
    data = pd.read_csv(annotated_data)
    functions(data)
    histogram(data)
    boxplot(data)
