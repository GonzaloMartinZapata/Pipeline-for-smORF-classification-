# Pipeline-for-smORF-classification

Genomes harbor an enormous number of Open Reading Frames (ORFs); however, determining which ones will be translated is not trivial, complicating the process of annotating protein-coding sequences. This challenge becomes more pronounced as the size of the ORFs decreases. The use of modern sequencing and proteomic technologies has allowed the observation that many organisms present small Open Reading Frames (smORFs) that are actively transcribed and translated. These smORFs have the potential to encode small-sized polypeptides, known as SEPs (SmORFs Encoded Polypeptides). 

It has been demonstrated that SEPs are involved in a wide range of processes in bacteria.  The exponential increase in the number of genomic sequences deposited in databases, coupled with the development of new machine learning methods, has led to an improvement in the bioinformatic prediction of smORFs that have the potential to be translated into SEPs. In fact, many genomes that have been recently updated in terms of annotation exhibit a variable number of annotated smORFs. However there is no much information on these genes in a large group of bacterias such as Plant Growth Promoting Rhizobacteria (PGPR), most of the bacterial plant pathogens, and certain genera responsible for causing diverse diseases in animals.

In this work, we conduct a comprehensive analysis of the current information already annotated about smORFs in the genomes of diverse bacterial species known to interact with eukaryotic hosts, focusing in particular on the analysis of smORFs from rhizobia, PGPRs and certain pathogens of plants or animals. These pipeline consists of the following steps:

![Schematic representation of the pipeline developed for SEP characterization.](https://raw.githubusercontent.com/GonzaloMartinZapata/Pipeline-for-smORF-classification-/main/Fig1.png)

**Figure 1:** Schematic representation of the pipeline developed for SEP characterization.

## 1- Download genomes for the genera/species of interest using NCBI’s command-line program datasets

The genomic sequences used in this work were downloaded for each taxon of interest using NCBI’s datasets program with the following command:

datasets download genome taxon "Bacteria of interest" --assembly-level complete --assembly-source RefSeq --exclude-atypical --include gbff --dehydrated --filename Bacteria_of_interest_dataset.zip

As a result, for each of the analyzed species, a large data package was downloaded as a dehydrated zip archive containing only metadata and the location of the data on NCBI servers. Once downloaded, the zip file was unzipped into a new folder and then “rehydrated” with the following command:

datasets rehydrate --directory “bacteria_of_interest”/

Next, we used a bash script ([rename.sh](rename.sh.txt)) provided by the NCBI to rename each file using the assembly accession number of the genome considered.

When more than 100 sequences met the selection criteria, the following analysis was performed. First, all the selected genomes were downloaded, and the FastANI software was used to calculate all the pairwise average nucleotide identity (ANI) values between a genome used as a reference and all the downloaded genomes ([FastANI tree script](Get_genomes_and_make_fastani_tree.ipynb)). Approximately fifty genomes were then selected by sampling the complete ANI value range.

Three particular cases that presented a larger number of genomic sequences were also considered: - **Bacillus**: We applied the same selection method as with the other genera but considered only genomes available from January 1, 2022, onward. This approach ensured a wide variety of species within the genus while maintaining the necessary variation. - **Pseudomonas** and **Klebsiella**: Most genomes were for a particular species. In these cases, approximately 25 genomes from the three most studied species were used. This approach was considered to provide results that more broadly represented each genus.


## 2- Extraction of the annotated SEP sequences and their associated information.

Two different python scripts were developed to work with the downloaded genomes. These scripts use Biopython {cita} to access the previously downloaded genbank files and extract the SEPs and the annotation information associated with them. The first one, that we called “annotated_smorfs” (https://github.com/GonzaloMartinZapata/Pipeline-for-smORF-classification-/blob/main/annotated_smorfs.py), extracts all the features that are translated and have a size of less than or equal to 70 amino acids. This script produces a fasta file foreach genome with all the coding sequences that meet our criteria. These files were used later to study the conservation of each smORF in the pangenome. The input for this script was a folder containing all the Genbank files of the considered species, and was runned with the following command

python annotated_smorfs.py -f “Bacteria_of_interest”

The second one, “annotated_info” (https://github.com/GonzaloMartinZapata/Pipeline-for-smORF-classification-/blob/main/annotated_info.py), works in a similar way. With the same input, this script generates a .csv file with the following information out of each smORF that was found: the accession number of the genome from which the smORF was extracted, the locus tag, the sequence, the length of that sequence, the product encoded for that smORF and, if it was available, the function of that SEP. Also this script takes some of this data to make three different graphics: a histogram illustrating the size distribution of all the smORFs found, a pie chart displaying the distribution of annotated products, and a box plot depicting the average size distribution per genome.

python annotated_info.py -f “Bacteria_of_interest”

## 3- Classification of smORFs by their degree of evolutionary conservation.

We used three different programs originally developed for pangenome analysis to evaluate smORF conservation along the different genomes (Orthofinder {ref}, Panaroo {ref} and Ortho GNC {ref}). To improve on speed we ran these programs with smORFs/SEPs sequences only, discarding all the other protein encoding genes. The commands and running modes used for each program were: 

For orthofinder,  orthofinder  –f ./ -S blast -og

We used Blast instead of diamond since it was reported to be at least 1-2% more accurate, and may be also more adequate for short peptides. In addition, the analysis was stopped after inferring the Orthogroups (-og option) as we were not interested in the phylogenetic study (for paralog discrimination) of the considered species.

For panaroo, the inputs were files in Prokka GFF3 format.  As we used only the smORFs sequences, a GFF3 file with the annotation information for smORFs only had to be made (https://github.com/GonzaloMartinZapata/Pipeline-for-smORF-classification-/blob/main/gbToGFF.py). For that purpose we used a custom script that “filtered” the Genbank files by size and then using the GFF python module converted it into the required format. Then, Panaroo was run with the following commands:

mkdir Panaroo_results
panaroo -i *.gff -o results –clean-mode strict

Lastly, Ortho GNC was run using the extracted protein sequences and the recommended settings.  

Once we had the results, a new python script was made to process the results of the three pangenome programs, and generate a csv file with all the gathered information (https://github.com/GonzaloMartinZapata/Pipeline-for-smORF-classification-/blob/main/results.py). First, we considered as “Conserved-50” all the smORFs that were present in 50% or more of the studied genomes. Next, we evaluated the presence/absence matrix given by each program, and classified the smORFs according to their conservation percentage. The smORFs were thus classified in six categories: “Probably not conserved” if the smORF was found in only one genome”, “Cloud” (0-15% presence in all considered genomes), “Shell” (15-50%), “Soft Core” (50-95%), “Core” (95-99%), and “conserved” (100%). This classification was later added to the .csv file generated by the “annotated_info” script. To do so, “results.py” script takes as input the “Orthogroups.txt” and “Orthogroups.GeneCount.tsv” files generated by Orthofinder, the “gene_data.csv” and “gene_presence_absence_roary.csv” files produced by Panaroo and the “OutputV2” generated by OrthoGNC. Then, it makes a “consensus” classification with the following criteria that is later added to the generated csv file: 
* If all programs classified a smORF into the same category, then that category was added to the results column
* If two out of three programs classified the smORF in the same category but the other one did not assign it to any category, or assigned it to a “neighbor” category (a category immediately before or immediately after the considered category), the consensus value was added to the results
* If the three results were different, or two were equal but the third was not classified in a “neighbor” category, the category of the lowest rank was added to the results column.

Those smORFs that were not assigned to a category by at least two of the programs were considered “Probably not conserved genes” 

For this script, the command used was: 

python results.py -og Orthogroups.txt -ofgc Orthogroups.GeneCount.tsv -pgd gene_data.csv -ppa gene_presence_absence_roary.csv -ognc OutputV2.txt -csv Bacteria_of_interes.csv -n xx

where n indicates the number of genomes considered for that species. 

## 4- Analysis of the obtained results

Finally, a last script (https://github.com/GonzaloMartinZapata/Pipeline-for-smORF-classification-/blob/main/analysis.py) was used to evaluate the obtained results. This script gives the number of smORFs found in the considered genre/species, the average number of smORFs per genome, a pie chart with the distribution according to the conservation rate among other data of interest. The command for this script is:

python analysis.py -r bacteria_of interest.csv -n xx 

The results obtained by this pipeline for all considered genera and species can be found in https://drive.google.com/drive/folders/13I8wmWN56fgk3GREeObXKApagEoQqQg1?usp=drive_link

