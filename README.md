# Pipeline-for-smORF-classification

Genomes harbor an enormous number of Open Reading Frames (ORFs); however, determining which ones will be translated is not trivial, complicating the process of annotating protein-coding sequences. This challenge becomes more pronounced as the size of the ORFs decreases. The use of modern sequencing and proteomic technologies has made it possible to observethat many organisms present small Open Reading Frames (smORFs) that are actively transcribed and translated. These smORFs have the potential to encode small-sized polypeptides, known as SEPs (SmORFs Encoded Polypeptides), that are involved in a wide range of biological processes. 
The exponential increase in the number of genomic sequences deposited in databases, coupled with the development of new machine learning methods, has led to an improvement in the bioinformatic prediction of smORFs that have the potential to be translated into SEPs. In fact, many genomes that have been recently updated in terms of annotation exhibit a variable number of annotated smORFs. However there is not much information about these genes in a large group of bacterias such as Plant Growth Promoting Rhizobacteria (PGPR), most of the bacterial plant pathogens, and certain genera responsible for causing diverse diseases in animals.
We developed a pipeline (Figure 1) to conduct a comprehensive analysis of the current information already annotated about smORFs in the genomes of diverse bacterial species known to interact with eukaryotic hosts, focusing in particular on the analysis of smORFs from rhizobia, PGPRs and certain pathogens of plants or animals. The pipeline consists of the following steps:

![Schematic representation of the pipeline developed for SEP characterization.](https://raw.githubusercontent.com/GonzaloMartinZapata/Pipeline-for-smORF-classification-/main/Fig1.png)

**Figure 1:** Schematic representation of the pipeline developed for SEP characterization.

## 1- Download genomes for the genera/species of interest using NCBI’s command-line program datasets

The genomic sequences can be downloaded for each taxon of interest using NCBI’s datasets program with the following command:

*datasets download genome taxon "Bacteria of interest" --assembly-level complete --assembly-source RefSeq --exclude-atypical --include gbff --dehydrated --filename Bacteria_of_interest_dataset.zip*

As a result, for each of the analyzed species, a large data package is downloaded as a dehydrated zip archive containing only the metadata and location of the data on NCBI servers. Once downloaded, the zip file can be unzipped into a new folder and then “rehydrated” with the following command:

*datasets rehydrate --directory “bacteria_of_interest”/*

Next, the rename bash script ([rename.sh](rename.sh.txt)) provided by the NCBI can be used to rename each file using the assembly accession number of the genome considered.

In the associated paper when more than 100 sequences met our selection criteria, the following analysis was performed. First, all the selected genomes were downloaded, and the FastANI software was used to calculate all the pairwise average nucleotide identity (ANI) values between a genome used as a reference and all the downloaded genomes ([FastANI tree script](Get_genomes_and_make_fastani_tree.ipynb)). Approximately fifty genomes were then selected by sampling the complete ANI value range.

## 2- Extraction of the annotated SEP sequences and their associated information.

Two different python scripts were developed to work with the downloaded genomes. These scripts use Biopython to access the previously downloaded genbank files and extract the SEPs and the annotation information associated with them. The first one, that we called “annotated_smorfs” ([Annotated_smorfs.py](annotated_smorfs.py)), extracts all the features that are translated and have a size of less than or equal to 70 amino acids. This script produces a fasta file foreach genome with all the coding sequences that meet our criteria. These files were used later to study the conservation of each smORF in the pangenome. The input for this script is a folder containing all the Genbank files of the considered species, and can be run with the following command:

*python annotated_smorfs.py -f “Bacteria_of_interest”*

The second one, “annotated_info” ([Annotated_info.py](annotated_info.py)), works in a similar way. With the same input, this script generates a .csv file with the following information out of each smORF that was found: the accession number of the genome from which the smORF was extracted, the locus tag, the sequence, the length of that sequence, the product encoded for that smORF and, if it was available, the function of that SEP (as annotated in the genbank file feature ‘function’ feature qualifier).  Also this script makes three different graphics: a histogram and a box plot illustrating the size distribution of the smORFs, and a pie chart displaying the functional distribution of the annotated products.

*python annotated_info.py -f “Bacteria_of_interest”*

## 3- Classification of smORFs by their degree of evolutionary conservation.

We use three different programs originally developed for pangenome analysis to evaluate smORF conservation along the different genomes (Orthofinder, Panaroo and Ortho GNC). To improve speed this pipeline is run with smORFs/SEPs sequences only, discarding all the other protein encoding genes. The commands and running modes used for each program are:

For orthofinder,  *orthofinder  –f ./ -S blast -og*

We use BLAST instead of DIAMOND since it was reported to be at least 1-2% more accurate, and may be also more adequate for short peptides. In addition, the analysis is stopped after inferring the Orthogroups (-og option) as we are not interested in the phylogenetic study (for paralog discrimination) of the considered species.

For panaroo, the inputs are files in Prokka GFF3 format. As we used only the smORFs sequences, a GFF3 file with the annotation information for smORFs only has to be made. To that end the provided script ([gbToGFF.py](gbToGFF.py)) can be used. This script filters the Genbank files by size and then using the GFF python module converts it into the required format. 

Then, Panaroo is run with the following commands:

*mkdir Panaroo_results*

*panaroo -i *.gff -o results –clean-mode strict*

Lastly, Ortho GNC is run using the extracted protein sequences and the recommended settings.  

## 4- Analysis of the obtained results

The results of the three pangenome programs are processed to generate consensus classification based on the degree of evolutionary conservation using the ([results.py](results.py)) script. First, the results given by each program are converted into a presence/absence matrix and the smORFs are classified into six categories according to their degree of conservation throughout the genomes: “Probably not conserved” if the smORF are found in only one genome, “Cloud” (0-15% presence in all considered genomes), “Shell” (15-50%), “Soft Core” (50-95%), “Core” (95-99%), and “Conserved” (100%). Next, a consensus classification is made with the following criteria: a) if all programs classify a smORF into the same category, then that category is assigned; b) if two programs classify the smORF in the same category but the third one does not assign it to any category, or assigns it to a “neighbor” category (a category immediately before or immediately after the considered category), the consensus value is assigned; c) if the three results are different, or two are equal but the third is not classified in a “neighbor” category, the category of the lowest rank is assigned; and d) those smORFs that are not assigned to a category by at least two of the programs are considered as “Probably not conserved genes”. 

This script takes as input the “Orthogroups.txt” and “Orthogroups.GeneCount.tsv” files generated by Orthofinder, the “gene_data.csv” and “gene_presence_absence_roary.csv” files produced by Panaroo and the “OutputV2” generated by OrthoGNC. 

To run the script the following command can be used:

*python results.py -og Orthogroups.txt -ofgc Orthogroups.GeneCount.tsv -pgd gene_data.csv -ppa gene_presence_absence_roary.csv -ognc OutputV2.txt -csv Bacteria_of_interes.csv -n xx*

where n indicates the number of genomes considered for that species. 

Finally, the script ([analysis.py](analysis.py)) is used to evaluate the obtained results. This script gives the number of smORFs found in the considered genera/species, the average number of smORFs per genome, and a pie chart with the distribution according to the conservation rate, among other data of interest. The command to run this script is:


*python analysis.py -r bacteria_of interest.csv -n xx*


The obtained tables for all genera and species considered in this work can be found in the [`Results_from_pipeline`](Results_from_pipeline/) directory.

The same tables allong with the graphs generated allong the pipeline are available at: 

https://drive.google.com/drive/folders/13I8wmWN56fgk3GREeObXKApagEoQqQg1?usp=drive_link

