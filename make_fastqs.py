#!/usr/bin/env python3
# Import csv and other modules
import argparse
import os
import glob
import sys
import shutil

script_usage = """A script for renaming entries in a multi-fasta file.
This script assumes sample ids in the fasta file are the string between ">" and the first "_" characters.
"""

parser = argparse.ArgumentParser(usage=script_usage)

parser.add_argument('-f', '--fastq', default='fastq_pass', help="directory containing fastq files")
parser.add_argument('-k', '--key', default="key_"+os.path.basename(os.getcwd())+".txt", help="Tab delimited file, where column 1 is the existing header and column 2 is the GISAID header, and column 3 is the Genbank header.")
parser.add_argument('-o', '--output', default="renamed_fastq", help="Directory of fastq files with new names")

args = parser.parse_args()

id_counter = 0
run_id=os.path.basename(os.getcwd())

# Check that input files exist
if not (os.path.exists(args.key)):
    sys.exit("\nERROR: An key file named %s doesnt exist\n" % (args.key))

if not os.path.isdir(args.fastq):
    sys.exit("\nERROR: A fastq dir named %s doesnt exist\n" % (args.fastq))

#Check for existing output
if os.path.isdir(args.output):
    sys.exit("\nERROR: A output dir named %s already exists\n" % (args.output))
else:
    os.mkdir(args.output)
    
rename_log = open(args.output+"/run_id_rename_log.tsv",'w')

# read in input file csv.DictReader
with open(args.key,'r') as id_file:
    for line in id_file:
        if line.count('\t') > 2:
            sys.exit("\nERROR: key file appears to contain more than 3 columns (%d tabs found) in line below:\n%s\n" % (line.count('\t'),line))
        #print(line)
        accession=line.split("\t")[0]
        seq_id=line.replace(' ','\t').split("\t")[2]
        fastq_path=args.fastq+"/"+accession+"_"+run_id+"*"
        source_list=glob.glob(fastq_path)
        dest_path=args.output+"/"+seq_id+".fastq.gz"
        
        #print(accession)
        print(seq_id)
        #print(fastq_path)
        
        if len(source_list) == 1:
            
            print("found it: %s" % source_list[0])
            shutil.copy(source_list[0],dest_path)
            rename_log.write(source_list[0]+"\t"+dest_path+"\n")
            id_counter=id_counter+1
        else:
            print("uh oh")
            print(glob.glob(fastq_path))
            sys.exit("found %s fastq files matching %s" % (len(source_list),fastq_path))

rename_log.close()

os.system('chmod -R 770 '+args.output)

print("Copied %d files\n" % id_counter)
