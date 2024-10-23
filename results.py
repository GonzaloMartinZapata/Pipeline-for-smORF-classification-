import pandas as pd
import os 
import csv
import argparse
parser= argparse.ArgumentParser (description= "takes the results files to classify each smORF")
parser.add_argument ("--orthogroups", "-og" ,dest= "orthogroups", required= True, help= "path to the file with the Orthofinder results" , nargs = 1) 
parser.add_argument ("--orthofinder_genecount", "-ofgc", dest= "orthofinder_genecount", required=True, help= "path to the genecount file of Orthofinder", nargs=1)
parser.add_argument ("--panaroo_gene_data", "-pgd", dest= "panaroo_gene_data" ,required=True, help= "gene_data file from Panaroo", nargs=1)
parser.add_argument ("--panaro_presence_absence", "-ppa", dest="panaroo_presence_absence", required=True, help= "gene_presence_absence_roary file of Panaroo", nargs=1)
parser.add_argument ("--ognc", "-ognc", dest="ognc", required=True, help="file with the OrtoGNC results", nargs=1)
parser.add_argument ("--bacteria_csvfile", "-csv", dest= "bacteria_csvfile", required=True, help= "path to the file with all the smORFs of a Bachteria", nargs=1)
parser.add_argument ("--number_of_genomes", "-n", dest="number_of_genomes", required=True, type=int, help="number of genomes where the smORFs were extracted", nargs=1 )
args= parser.parse_args ()

csv_file= args.bacteria_csvfile [0]

#   Fist step: add the Orthofinder results to the CSV
print ("Processing Orthofinder Results")
orthogroups= args.orthogroups [0]
orthologs_data = []

with open (orthogroups, "r") as in_handle: 
    for line in in_handle:
        division= line.strip().split(":")
    
        group = division[0].strip()
        genes = division[1].strip().split(' ')
        dict_orthogroups = {"Orthogroup": group}
        for i, gen in enumerate(genes, start=1):
            dict_orthogroups[f"Gene_{i}"] = gen

        orthologs_data.append(dict_orthogroups)

df_orthogroups = pd.DataFrame(orthologs_data)

#Classify the different Orthogroups in the gene_count file
genecount = args.orthofinder_genecount[0]

df_gc= pd.read_csv (genecount, sep="\t")
df_gc= df_gc.drop (["Total"], axis=1)

for column in df_gc.columns [1:]:
    df_gc [column]= df_gc[column].apply(lambda x: True if x > 0 else False)

true_count_per_row = df_gc.iloc[:, 1:].sum(axis=1)

df_gc["presence_percentage"] = ((true_count_per_row/len (df_gc.columns[1:])) *100) 

bins = [0, 15, 50, 95, 99, 100]
labels = ["Cloud", "Shell", "Soft core", "Core", "Conserved"]

df_gc ["Classification_of"] = pd.cut(df_gc["presence_percentage"], bins=bins, labels=labels)

df_orthofinder = pd.merge(df_orthogroups, df_gc[["Orthogroup", "Classification_of"]], on="Orthogroup", how="left")


bacteria_df = pd.read_csv(csv_file, sep=",", quotechar='"', quoting=csv.QUOTE_MINIMAL)

bacteria_df ["Orthogroup"]= ""
bacteria_df ["Classification_of"]= ""

for index, row in bacteria_df.iterrows():
    locus_tag_value = row["locus_tag"]
    matching_row = df_orthofinder[df_orthofinder.isin([locus_tag_value]).any(axis=1)]
    if not matching_row.empty:
        bacteria_df.at[index, "Orthogroup"] = matching_row["Orthogroup"].values[0]
        bacteria_df.at[index, "Classification_of"] = matching_row["Classification_of"].values[0]

print ("Orthofinder results added")

#   Part2: Add the Panaroo results

print ("Processing Panaroo Results")
gene_data= args.panaroo_gene_data [0]
presence_absence= args.panaroo_presence_absence [0]

df_gd = pd.read_csv(gene_data)
df_pa= pd.read_csv(presence_absence)

no_isolates = df_pa["No. isolates"]
total_of_samples = args.number_of_genomes [0]

df_pa["presence_percentage"] = (no_isolates / total_of_samples) * 100
bins = [0, 15, 50, 95, 99, 100]
labels = ["Cloud", "Shell", "Soft core", "Core", "Conserved"]

df_pa["Classification_panaroo"] = pd.cut(df_pa["presence_percentage"], bins=bins, labels=labels)

gene_values_list = []
classification_values_list = []

for index, row in df_gd.iterrows():
    annotation_id_value = row["annotation_id"]
    gff_files_value = row["gff_file"]

    if gff_files_value in df_pa.columns:
        gene_value = df_pa.loc[df_pa[gff_files_value] == annotation_id_value, "Gene"].values
        classification_value = df_pa.loc[df_pa[gff_files_value] == annotation_id_value, "Classification_panaroo"].values

        gene_values_list.append(gene_value[0] if len(gene_value) > 0 else None)
        classification_values_list.append(classification_value[0] if len(classification_value) > 0 else None)
    else:
        gene_values_list.append(None)
        classification_values_list.append(None)

df_gd["Gene"] = gene_values_list
df_gd["Clasification_panaroo"] = classification_values_list

gene_values_list = []
clasification_panaroo_values_list = []

for index, row in bacteria_df.iterrows():
    genome_value = row["genome"]
    prot_seq_value = row["prot_seq"]

    matching_rows = df_gd[(df_gd["gff_file"] == genome_value) & (df_gd["prot_sequence"] == prot_seq_value)]

    if len(matching_rows) > 0 and all(matching_rows["gff_file"] == genome_value) and all(matching_rows["prot_sequence"] == prot_seq_value):
        gene_values = matching_rows["Gene"].tolist()
        clasification_panaroo_values = matching_rows["Clasification_panaroo"].tolist()

        gene_str = ', '.join(filter(None, gene_values))
        clasification_panaroo_str = ', '.join(filter(None, clasification_panaroo_values))

        bacteria_df.at[index, "Gene"] = gene_str
        bacteria_df.at[index, "Classification_panaroo"] = clasification_panaroo_str
    else:
        bacteria_df.at[index, "Gene"] = None
        bacteria_df.at[index, "Classification_panaroo"] = None

print ("Panaroo results added")

#   part3: Add the OGNC results

print ("Processing OGNC Results")

ognc= args.ognc [0]

with open (ognc, "r") as in_handle:
    lines= in_handle.readlines()

data= []

for line in lines:
    parts = line.strip().split('\t', 1)
    # Extract the first locus_tag in the line as it is the considered gene
    first_locus_tag = parts[0]
    
    # extract the other locus tag coma separated in the second part. These are the orthologs 
    if len(parts) > 1:
        additional_locus_tags = parts[1].split(',')
        
        # group the locus tags that are from the same genome (they are equal in the part before the "_")
        grouped_genes = {}
        for locus_tag in additional_locus_tags:
            genome = locus_tag.split('_')[0]
            if genome not in grouped_genes:
                grouped_genes[genome] = [locus_tag]
            else:
                grouped_genes[genome].append(locus_tag)
        
        row_data = [first_locus_tag]
        for genome_genes in grouped_genes.values():
            row_data.append(','.join(genome_genes))
        
        data.append(row_data)
    else:
        data.append([first_locus_tag])

df_ognc = pd.DataFrame(data)


number_of_genomes= args.number_of_genomes [0]

max_columns = max(len(row) for row in data)

df_ognc["presence_percentage"] = df_ognc.apply(lambda row: row.count() / max_columns * 100, axis=1)

bins = [0, 15,50, 95, 99, 100]
labels = ["Cloud", "Shell", "Soft core", "Core", "Conserved"]

df_ognc["Classification_ognc"] = pd.cut(df_ognc["presence_percentage"], bins=bins, labels=labels)

df_ognc = df_ognc.rename(columns={df_ognc.columns[0]: "locus_tag"})

bacteria_df = pd.merge(bacteria_df, df_ognc[["locus_tag", "Classification_ognc"]], how="left", left_on="locus_tag", right_on="locus_tag")

print ("OGNC results added")

#   Part 4, make a "consensus" result
print ("Processing results")

priority_order = ["Cloud", "Shell", "Soft core", "Core", "Conserved"]

df= bacteria_df
# Handle cases where "Classification_panaroo" has two values separated by a comma
df["Classification_panaroo"] = df["Classification_panaroo"].apply(lambda x: min(x.split(','), key=lambda val: priority_order.index(val.strip())) if pd.notna(x) and ',' in str(x) else x)

df["Results_column"] = None  

#Apply different conditions to classify smORFs

# Conditions 1, 2 and 3: Only one program classified the given smORF thus is considered as "probably not conserved"

condition_1 = df["Classification_of"].notna() & df["Classification_panaroo"].isna() & df["Classification_ognc"].isna()
df.loc[condition_1, "Results_column"] = "Probably not conserved"

condition_2 = df["Classification_of"].isna() & df["Classification_panaroo"].notna() & df["Classification_ognc"].isna()
df.loc[condition_2, "Results_column"] = "Probably not conserved"

condition_3 = df["Classification_of"].isna() & df["Classification_panaroo"].isna() & df["Classification_ognc"].notna()
df.loc[condition_3, "Results_column"] = "Probably not conserved"

# Apply further conditions in order, only if Results_column is None

# Condition 4: All three programs classified in the same category
condition_4 = (df["Classification_of"] == df["Classification_panaroo"]) & (df["Classification_of"] == df["Classification_ognc"]) & df["Results_column"].isna()
df.loc[condition_4, "Results_column"] = df["Classification_of"]

# Condition 5, 6 and 7: Two programs classified the same, and the third did not assign a category
condition_5 = (df["Classification_of"] == df["Classification_panaroo"]) & df["Classification_ognc"].isna() & df["Results_column"].isna()
df.loc[condition_5, "Results_column"] = df["Classification_of"]

condition_6 = (df["Classification_of"] == df["Classification_ognc"]) & df["Classification_panaroo"].isna() & df["Results_column"].isna()
df.loc[condition_6, "Results_column"] = df["Classification_of"]

condition_7 = (df["Classification_panaroo"] == df["Classification_ognc"]) & df["Classification_of"].isna() & df["Results_column"].isna()
df.loc[condition_7, "Results_column"] = df["Classification_panaroo"]

# Condition 8: Two classifications are equal and the third is a neighbor
neighbours = {"Cloud": ["Shell"],
              "Shell": ["Cloud", "Soft core"],
              "Soft core": ["Shell", "Core"],
              "Core": ["Soft core", "Conserved"],
              "Conserved": ["Core"]}

for index, row in df.iterrows():
    if pd.isna(row["Results_column"]):
        if row["Classification_of"] == row["Classification_panaroo"] and row["Classification_of"] != row["Classification_ognc"]:
            if row["Classification_ognc"] in neighbours.get(row["Classification_of"], []):
                df.at[index, "Results_column"] = row["Classification_of"]
        elif row["Classification_of"] == row["Classification_ognc"] and row["Classification_of"] != row["Classification_panaroo"]:
            if row["Classification_panaroo"] in neighbours.get(row["Classification_of"], []):
                df.at[index, "Results_column"] = row["Classification_of"]
        elif row["Classification_panaroo"] == row["Classification_ognc"] and row["Classification_panaroo"] != row["Classification_of"]:
            if row["Classification_of"] in neighbours.get(row["Classification_panaroo"], []):
                df.at[index, "Results_column"] = row["Classification_panaroo"]

# Condition 9: All three results are different; take the lowest value
condition_9 = (df["Classification_of"] != df["Classification_panaroo"]) & (df["Classification_of"] != df["Classification_ognc"]) & (df["Classification_panaroo"] != df["Classification_ognc"]) & df["Results_column"].isna()
for index in df.loc[condition_9].index:
    row = df.loc[index]
    cat_row = pd.Categorical([row["Classification_of"], row["Classification_panaroo"], row["Classification_ognc"]], categories=priority_order, ordered=True)
    df.at[index, "Results_column"] = cat_row.min()

# Condition 10: Two classifications are equal and the third is not a neighbor; take the lowest value
condition_10 = ((df["Classification_of"] == df["Classification_panaroo"]) & (df["Classification_of"] != df["Classification_ognc"]) | 
               (df["Classification_of"] == df["Classification_ognc"]) & (df["Classification_of"] != df["Classification_panaroo"]) |
               (df["Classification_panaroo"] == df["Classification_ognc"]) & (df["Classification_panaroo"] != df["Classification_of"])) & df["Results_column"].isna()
for index in df.loc[condition_10].index:
    row = df.loc[index]
    cat_row = pd.Categorical([row["Classification_of"], row["Classification_panaroo"], row["Classification_ognc"]], categories=priority_order, ordered=True)
    df.at[index, "Results_column"] = cat_row.min()

# Condition 11: None of the programs assigned a category
condition_11 = df["Classification_of"].isna() & df["Classification_panaroo"].isna() & df["Classification_ognc"].isna() & df["Results_column"].isna()
df.loc[condition_11, "Results_column"] = "Probably not conserved"

results_directory= "results"
if not os.path.exists (results_directory):
    os.mkdir (results_directory)

df.to_csv (os.path.join (results_directory, csv_file), index=False)
