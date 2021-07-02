#!/usr/bin/env python

# std libraries
import os
import sys
import time
import csv
import subprocess
import argparse
from Bio import SeqIO
import pandas as pd


class Genome:
    def __init__(self):
        pass

class Gene:
    def __init__(self):
        pass
    
class Reference:
    def __init__(self):
        pass

class Target:
    def __init__(self):
        pass


# creating the parser
parser = argparse.ArgumentParser()

# parser arguments
parser.add_argument('--fna', type = str, required = True, help = 'Input fasta file or path to the file')
parser.add_argument('--db', type = str, required = True, help = 'Database file or path to database')
parser.add_argument('--threads', '-t', type = int, required = True, help = 'Number of CPU threads; used for diamond.')
# parser.add_argument('--diamond_out', type = str, required = True, help = 'Diamond output file name')

args = parser.parse_args()


# get basenames for input files
inputbase = os.path.basename(f'{args.fna}') # currently have set to zero since the directory contains two files
inputbase = os.path.splitext(inputbase)[0] # the output is a tuple that where the first value the name   

db_base = os.path.basename(f'{args.db}')
db_base = os.path.splitext(db_base)[0]



########################## gene predictions with Prodigal ########################## 

def run_prodigal(input):
    """
    Running gene predictions with Prodigal
    1. Using anonymous mode (-p meta has been depreciated)
    2. Output as genbank format and protein sequences (faa)
    """
    
    # creating the shell script command
    prod_cmd = "prodigal "
    prod_cmd += f"-i {input} " # input is fna 
    prod_cmd += f"-o {inputbase}_prodigal.gbk " 
    prod_cmd += f"-a {inputbase}_prodigal.faa " # this is an output file used
    prod_cmd += "-p meta"
    
    # running the command
    p = subprocess.call(prod_cmd, shell = True) # for testing printing the command


########################## Make diamond database ########################## 

def make_diamonddb(database):
    """
    Creating a database for diamond
    """

    db_base = os.path.basename(f'{database}')
    db_base = os.path.splitext(db_base)[0]
    # creating the shell script command
    diamonddb_cmd = 'diamond '
    diamonddb_cmd += 'makedb '
    diamonddb_cmd += f'--in {database} '
    diamonddb_cmd += f'-d {db_base}'
    
    # running the command
    p = subprocess.call(diamonddb_cmd, shell = True) # for testing printing the command


########################## Running diamond ########################## 

def run_diamond(input_faa, diamond_db, output):
    """
    Running diamond
    """
    
    # creating the shell script with the command
    diamond_cmd = "diamond "
    diamond_cmd += "blastp "
    diamond_cmd += f"-q {input_faa} " # use output from prodigal as input for diamond
    diamond_cmd += f"-d {diamond_db} "
    diamond_cmd += "-e 0.0001 " # maximum e-value of 10^-4
    diamond_cmd += "--id 35 " # minimum identity% coverage greater than or equal to 35%
    diamond_cmd += "--query-cover 50 " # query coverage greater than or equal to 50%
    diamond_cmd += "--subject-cover 50 " # reference coverage greater than or equal to 50%
    diamond_cmd += "-k 0 " # max number of target seq for align (default = 25)
    diamond_cmd += f"-p {args.threads} "
    diamond_cmd += f"-o {output}_dmnd_out.tsv"
    
    # running the command
    p = subprocess.call(diamond_cmd, shell = True) # for testing printing the command

########################## Calculating PPQ/gPPQ Score ########################## 


def calc_PPQ(virus_hits, plasmid_hits, totalvirus, totalplasmid):
    """
    Function calculates the PPQ scores from Pfeifer et.al. 2021.
    """
    try:
        return (int(virus_hits) / int(totalvirus)) / ((int(virus_hits) / int(totalvirus))+(int(plasmid_hits) / int(totalplasmid)))
    except ZeroDivisionError:
        return None

    
def calc_gPPQ(sum_ppq, num_genes_wppq):
    """
    Function calculates the gPPQ scores from Pfeifer et.al. 2021.
    """
    try:
        return ((sum_ppq) / ((num_genes_wppq)))
    except ZeroDivisionError:
        return None


def main():
    
    startTime = time.time()
    run_prodigal(args.fna)
    make_diamonddb(args.db)
    run_diamond(input_faa=f'{inputbase}_prodigal.faa', diamond_db=f'{db_base}.dmnd', output=f'{inputbase}')
    
    ########################## genome object ##########################

    ## collect genome information
    genomes = {}
    for header in SeqIO.parse(args.fna, 'fasta'):
        genome = Genome()
        genome.id = header.id.split()[0]
        genome.len = len(str(header.seq).upper())
        genome.genes = 0
        genome.num_genes_ppq = 0
        genome.ppq_sum = None
        genome.gppq = None
        genomes[genome.id] = genome

    ########################## genes object ##########################
    
    genes = {}

    for header in SeqIO.parse(f'{inputbase}_prodigal.faa', 'fasta'):
        gene = Gene()
        gene.id = header.id.split()[0]
        gene.genome_id = header.id.rsplit('_',1)[0]
        gene.len = len(str(header.seq).upper())
        gene.virus = set([])
        gene.plasmid = set([])
        gene.ppq = None
        genes[gene.id] = gene
        genomes[gene.genome_id].genes += 1


    ########################## reference ##########################
    
    reference = {}
    total_dbtype = []
    reader = csv.DictReader(open('mergedDB_metadata.csv')) # add arg here
    
    for row in reader:
        ref = Reference()
        ref.gene_id = row['sequence_header']
        ref.database = row['database']
        ref.genome_id = row['genome_id']
        reference[ref.gene_id] = ref
        total_dbtype.append(row['database'])

    totalvirus = total_dbtype.count('virus')
    totalplasmid = total_dbtype.count('plasmid')
    
    del total_dbtype
    del reader
    print(f'Total virus database count: {totalvirus} \nTotal plasmid database count: {totalplasmid}')


    ########################## target ##########################
    
    dmnd_headers = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]


    # pythonic way of adding header to a tsv
    # https://stackoverflow.com/a/50129816
    with open(f'{inputbase}_dmnd_out.tsv', newline='') as f_input, open(f'{inputbase}_dmnd_out.csv', 'w', newline='') as f_output:
        r = csv.reader(f_input, delimiter='\t')
        w = csv.writer(f_output, delimiter=',') # convert to csv 

        w.writerow(dmnd_headers)
        w.writerows(r)


    targets = {}
    reader =  csv.DictReader(open(f'{inputbase}_dmnd_out.csv'))
    for row in reader: 
        target = Target()
        target.gene_id = row['qseqid']
        target.query_genome_id = reference[target.gene_id].genome_id
        target.match_id = row['sseqid']
        target.match_genome_id = reference[target.match_id].genome_id
        target.type = reference[target.match_id].database
        targets[target.gene_id] = target
    
    ##### Counting hits 
    for trgt in targets.keys():
        if targets[trgt].query_genome_id == targets[trgt].match_genome_id:
            continue
        if targets[trgt].type == 'virus':
            genes[targets[trgt].gene_id].virus.add(targets[trgt].match_genome_id)
        elif targets[trgt].type == 'plasmid':
            genes[targets[trgt].gene_id].plasmid.add(targets[trgt].match_genome_id)
        else:
            None


    ########################## PPQ/gPPQ ##########################

    for genome_id in genomes.items():
        if genomes[genome_id[0]].genes >= 10:
            g_id = genomes[genome_id[0]].id
            ppq_scores = []
            
            for gene_id in genes.items():
                if genes[gene_id[0]].genome_id == g_id:
                    vhits = genes[gene_id[0]].virus
                    phits = genes[gene_id[0]].plasmid
                    ppq = calc_PPQ(virus_hits=vhits,plasmid_hits=phits, totalvirus=totalvirus, totalplasmid=totalplasmid)
                    genes[gene_id[0]].ppq = ppq
                    if ppq is not None:
                        genomes[genome_id[0]].num_genes_ppq += 1
                        ppq_scores.append(ppq)
                    
            genomes[genome_id[0]].ppq_sum = sum(ppq_scores)
        if genomes[genome_id[0]].ppq_sum is not None:
            ppq_s = genomes[genome_id[0]].ppq_sum
            num_ppqgenes = genomes[genome_id[0]].num_genes_ppq
            gppq = calc_gPPQ(sum_ppq = ppq_s, num_genes_wppq = num_ppqgenes)
            genomes[genome_id[0]].gppq = gppq

    ### writing gPPQ scores to a csv file ###
    id_list = []
    num_genes = []
    num_ppq_genes = []
    genome_len = []
    gPPQ_scores = []
    for genome_id in genomes.items():
        id_list.append(genomes[genome_id[0]].id)
        num_genes.append(genomes[genome_id[0]].genes)
        num_ppq_genes.append(genomes[genome_id[0]].num_genes_ppq)
        genome_len.append(genomes[genome_id[0]].len)
        gPPQ_scores.append(genomes[genome_id[0]].gppq)

    out_df = pd.DataFrame({'genome_id':id_list, 
                           'number_genes':num_genes,
                           'num_genes_w_ppq':num_ppq_genes,
                           'genome_length':genome_len,
                           'gPPQ':gPPQ_scores})
    
    out_df.to_csv(f'{inputbase}_gPPQ_scores.csv', index=False)
    runtime = round((time.time() - startTime),2)
    print(f"Run time: {runtime} seconds")


        
if __name__ == "__main__":
    main()


