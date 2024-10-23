#!/usr/bin/env python3
'''
Converts Gb to GFF3 using Biopython and BCBio.
Extract Gb features to GFF3
Requires: BCBio - pip install bcbio-gff
'''
import sys, glob, os, argparse
from BCBio import GFF
from Bio import SeqIO

#functions

def parseArgs():
    '''
    Argument parsing is done.
    '''
    parser = argparse.ArgumentParser(description='Converts gb files to GFF3. Extracts a list of features to GFF3.')
    parser.add_argument("-g", "--genomes",help="Genomes to convert in Genbank format. <Space separated list of genbank files. Also accepts *.gb >", dest="genomes", nargs='+', action='store', required=True)
    parser.add_argument("-e", "--extract",help="Set on for feature extraction. <boolean>", dest="extract", action='store_true', required=False)
    parser.add_argument("-s", "--size",help="Size cutoff <INT>", dest="size", action='store', required=False, default = 210)
    parser.add_argument("-l", "--locus_tags",help="A file containing a list of locus_tags to extract. <File, one locus tag on each line>", dest="lt_file", action='store', required=False)
    args = parser.parse_args()
    return args

def filter_gb(lt_list,gb,size,extract_from_list):
    output = f'{os.path.splitext(gb)[0]}_filtered{os.path.splitext(gb)[1]}'
    r_list = []
    if extract_from_list:
        for r in SeqIO.parse(gb,'gb'):
            feature_list = []
            for feature in r.features:
                if feature.type == 'CDS':
                    if ('locus_tag' in feature.qualifiers) and ('translation' in feature.qualifiers):
                        if feature.qualifiers['locus_tag'][0] in lt_list:
                            feature_list.append(feature)
            r.features = feature_list
            r_list.append(r)
        SeqIO.write(r_list,output,'gb')
    else:
        lt_list = []
        for r in SeqIO.parse(gb,'gb'):
            feature_list = []
            for feature in r.features:
                if feature.type == 'CDS':
                    if ('locus_tag' in feature.qualifiers) and ('translation' in feature.qualifiers):
                        ltag = feature.qualifiers['locus_tag'][0]
                        feature_size = int(feature.location.end-feature.location.start)
                        if feature_size <= size:
                            feature_list.append(feature)
                            lt_list.append(ltag)
            r.features = feature_list
            r_list.append(r)
        SeqIO.write(r_list,output,'gb')
    return(output,lt_list)

def convert_to_GFF(f,extract, lt_list,size,extract_from_list):
    nf = os.path.splitext(f)[0]+'.gff'
    with open(f, 'r') as fh, open(nf,'+w') as nfh:
        if extract:
            edited_f, lt_list = filter_gb(lt_list,f,size,extract_from_list)
            GFF.write(SeqIO.parse(edited_f, "genbank"), nfh)
        else:
            lt_list = []
            GFF.write(SeqIO.parse(f, "genbank"), nfh)
        nfh.write('##FASTA\n')
        for record in SeqIO.parse(fh, "genbank"):
            SeqIO.write(record, nfh, 'fasta')
    return(lt_list)

if __name__ == '__main__':
    args = parseArgs()
    gbfiles = args.genomes
    filesToConvert = gbfiles
    #
    if '-' in gbfiles:
        sys.exit(f"ERROR: Incorrect options '-'")
    if len(gbfiles) == 1:
        if '*.' in gbfiles[0]:
            filesToConvert=glob.glob(gbfiles[0])
    #
    extract = args.extract
    size = int(args.size)
    extract_from_list = False
    lt_list = []
    if extract:
        if (not len(gbfiles) == 1) and (args.lt_file):
            sys.exit('ERROR: Feature extraction from locus_tag list only works with one genome file...')
        if not args.lt_file:
            print(f'Extraction of features of less than {size} nucleotides.')
        elif not os.path.exists(args.lt_file):
            print(f'ERROR: File not found. Extraction of features of less than {size} nucleotides.')
        else:
            extract_from_list = True
            with open(args.lt_file,'r') as ltf:
                for line in ltf:
                    lt_list.append(line.strip())
    else:
        if args.lt_file:
            sys.exit('ERROR: Locus_tag file passed without the extract option. Please set -e option...')
    #
    # Extract selected features and convert to GFF3
    for f in filesToConvert:
        if os.path.exists(f):
            print(f'Processing {f}...')
            lt_list = convert_to_GFF(f,extract, lt_list,size,extract_from_list)
            if extract:
                print(f'A total of {len(lt_list)} smORFs were extracted...')
        else:
            print(f'{f} not found, skipping to next file...')