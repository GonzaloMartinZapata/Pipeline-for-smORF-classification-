# Pipeline-for-smORF-classification

Genomes harbor an enormous number of Open Reading Frames (ORFs); however, determining which ones will be translated is not trivial, complicating the process of annotating protein-coding sequences. This challenge becomes more pronounced as the size of the ORFs decreases. The use of modern sequencing and proteomic technologies has allowed the observation that many organisms present small Open Reading Frames (smORFs) that are actively transcribed and translated. These smORFs have the potential to encode small-sized polypeptides, known as SEPs (SmORFs Encoded Polypeptides). 

It has been demonstrated that SEPs are involved in a wide range of processes in bacteria. However there is no much information on these genes in a large group of bacterias such as Plant Growth Promoting Rhizobacteria (PGPR), most of the bacterial plant pathogens, and certain genera responsible for causing diverse diseases in animals.

## Methods

1-  Genomic sequences used in this study
the used genomes were downloaded using NCBI´s command-line program "datasets" using the following command:

datasets download genome taxon "Bacteria of interest" --assembly-level complete  --assembly-source RefSeq --exclude-atypical --include gbff --dehydrated --filename Bacteria_of_interest_dataset.zip

where “Bacteria of interest" corresponds to the different genus and species. As a result, for each of the analyzed  species a large data package was downloaded as a dehydrated zip archive that contained only metadata and the location of the data on NCBI servers. Once downloaded, the zip file was unzipped in a new folder and then “rehydrated” with the following command:

datasets rehydrate  --directory “bacteria_of_interest”/

Next, we used a bash script (poner nombre del script o algo que lo identifique) provided by the NCBI to rename each file as the assembly accession number of the genome considered. 

2- Extraction of the annotated SEP sequences and their associated information.

Two different python scripts were developed to work with the downloaded genomes. These scripts use Biopython {cita} to access the previously downloaded genbank files and extract the SEPs and the annotation information associated with them. The first one, that we called “annotated_smorfs”, extracts all the features that are translated and have a size of less than or equal to 70 amino acids. This script produces a fasta file foreach genome with all the coding sequences that meet our criteria. These files were used later to study the conservation of each smORF in the pangenome. The input for this script was a folder containing all the Genbank files of the considered species, and was runned with the following command
python annotated_smorfs.py -f “Bacteria_of_interest”

The second one, “annotated_info”, works in a similar way. With the same input, this script generates a .csv file with the following information out of each smORF that was found: the accession number of the genome from which the smORF was extracted, the locus tag, the sequence, the length of that sequence, the product encoded for that smORF and, if it was available, the function of that SEP. Also this script takes some of this data to make three different graphics: a histogram illustrating the size distribution of all the smORFs found, a pie chart displaying the distribution of annotated products, and a box plot depicting the average size distribution per genome.

python annotated_info.py -f “Bacteria_of_interest”

3- Classification of smORFs by their degree of evolutionary conservation.
We used three different programs originally developed for pangenome analysis to evaluate smORF conservation along the different genomes (Orthofinder {ref}, Panaroo {ref} and Ortho GNC {ref}). To improve on speed we ran these programs with smORFs/SEPs sequences only, discarding all the other protein encoding genes. The commands and running modes used for each program were: 
For orthofinder,  orthofinder  –f ./ -S blast -og
We used Blast instead of diamond since it was reported to be at least 1-2% more accurate –and may be also more adequate for short peptides–, although with a runtime of approximately 20x longer. The -og option stopped the analysis after inferring the Orthogroups as we were not interested in the phylogenetic study of the considered species. 
For panaroo, the inputs were files in Prokka GFF3 format.  As we used only the smORFs sequences, a GFF3 file with the annotation information for smORFs only had to be made. For that purpose we used a custom script (subir al github) that “filtered” the Genbank files by size and then using the GFF python module converted it into the required format. Then, Panaroo was run with the following commands:

mkdir Panaroo_results
panaroo -i *.gff -o results –clean-mode strict

Lastly, Ortho GNC was run using the extracted protein sequences and the recommended settings.  
Once we had the results, a new python script was made to process the results of the three pangenome programs, and generate a csv file with all the gathered information. First, we considered as “probably conserved” all the smORFs that were present in 50% or more of the studied genomes. Next, we evaluated the presence/absence matrix given by each program, and classified the smORFs according to their conservation percentage. The smORFs were thus classified in six categories: “Probably not conserved” if the smORF was found in only one genome”, “Cloud” (0-15% presence in all considered genomes), “Shell” (15-50%), “Soft Core” (50-95%), “Core” (95-99%), and “conserved” (100%). This classification was later added to the .csv file generated by the “annotated_info” script. To do so, “results.py” script takes as input the “Orthogroups.txt” and “Orthogroups.GeneCount.tsv” files generated by Orthofinder, the “gene_data.csv” and “gene_presence_absence_roary.csv” files produced by Panaroo and the “OutputV2” generated by OrthoGNC. Then, it makes a “consensus” classification with the following criteria that is later added to the generated csv file: 
If all programs classified a smORF into the same category, then that category was added to the results column
If two out of three programs classified the smORF in the same category but the other one did not assign it to any category, or assigned it to a “neighbor” category (a category immediately before or immediately after the considered category), the consensus value was added to the results
If the three results were different, or two were equal but the third was not classified in a “neighbor” category, the category of the lowest rank was added to the results column.
Those smORFs that were not assigned to a category by at least two of the programs were considered “Probably not conserved genes” 

For this script, the command used was: 

python results.py -og Orthogroups.txt -ofgc Orthogroups.GeneCount.tsv -pgd gene_data.csv -ppa gene_presence_absence_roary.csv -ognc OutputV2.txt -csv Bacteria_of_interes.csv -n xx

where n indicates the number of genomes considered for that species. 

4-Bioinformatic characterization of hypothetical SEPs
We aimed at the functional characterization of SEPs that were classified as “probably conserved” and also were annotated as hypothetical proteins in selected strains of interest. For that purpose -select_genome.py- was developed. This script creates a fasta file with all the SEPs that meet this criteria for a chosen organism. The command-line used was:
python select_genome.py -csv “csv file of the considered organism” -g accession number of the strain to be studied. 
This fasta file was used later for characterization using the following programs:
Pannzer2 web server was used for the annotation of GO terms (the batch processing was used considering that we were using more than 10 sequences per study)
Batch web cd search tool was used to search for conserved domains, the used databases were CDD, COG, PFAM and Tigr. The expected value threshold was of 0.01
TMHMM-2.0 was used to search for transmembrane regions. 
SignalP-6.0 and Phobius were used to search for signal peptides.
IUPRED3 was used to search of disordered regions 
s4pred was used for the prediction of secondary structure. After installing the program, the command used was: python run_model.py  <fasta file with hypothetical proteins> > Results.ss2  