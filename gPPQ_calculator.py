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
import numpy as np


class Genome:
    def __init__(self):
        self.genes = []
        self.num_genes_ppq = 0
        self.ppq_sum = None
        self.gppq = None
        
    def calc_gppq(self, genes):
        """
        Function calculates the gPPQ scores from Pfeifer et.al. 2021.
        """
        self.ppq_scores = [genes[gene_id].ppq for gene_id in self.genes if genes[gene_id].ppq is not None]
        self.num_genes = len(self.genes)
        self.num_ppq = len(self.ppq_scores)
        self.gppq = np.mean(self.ppq_scores) if len(self.ppq_scores) > 0 else None
        
class Gene:
    def __init__(self):
        self.virus_hits = set([])
        self.plasmid_hits = set([])
        
    def calc_ppq(self, totalvirus, totalplasmid):
        """
        Function calculates the PPQ scores from Pfeifer et.al. 2021.
        """
        virus_hit_rate = 1.0 * len(self.virus_hits) / totalvirus
        plasmid_hit_rate = 1.0 * len(self.plasmid_hits) / totalplasmid
        if virus_hit_rate + plasmid_hit_rate == 0:
            self.ppq = None
        else:
            self.ppq = virus_hit_rate / (virus_hit_rate + plasmid_hit_rate)
    
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
parser.add_argument('--db', type = str, required = True, help = 'Path to database directory containing: proteins.dmnd,  proteins.faa,  proteins.tsv')
parser.add_argument('--threads', type = int, default = 1, help = 'Number of CPU threads; used for diamond')
parser.add_argument('--out', type = str, required = True, help = 'Output directory')
args = parser.parse_args()

# get basenames for input files
if not os.path.exists(args.out): os.makedirs(args.out)
inputbase = os.path.basename(f'{args.fna}') # currently have set to zero since the directory contains two files
inputbase = os.path.splitext(inputbase)[0] # the output is a tuple that where the first value the name   

db_base = os.path.basename(f'{args.db}')
db_base = os.path.splitext(db_base)[0]



########################## gene predictions with Prodigal ########################## 

def run_prodigal(fna, out):
    """
    Running gene predictions with Prodigal
    1. Using anonymous mode (-p meta has been depreciated)
    2. Output as genbank format and protein sequences (faa)
    """
    
    # creating the shell script command
    prod_cmd = "prodigal "
    prod_cmd += f"-i {fna} " # input is fna
    prod_cmd += f"-a {out}/proteins.faa " # this is an output file used
    prod_cmd += "-p meta"
    prod_cmd += "> /dev/null"
    
    # running the command
    p = subprocess.call(prod_cmd, shell = True) # for testing printing the command

########################## Running diamond ########################## 

def run_diamond(faa, db, out, threads):
    """
    Running diamond
    """
    
    # creating the shell script with the command
    diamond_cmd = "diamond blastp "
    diamond_cmd += f"-q {faa} " # use output from prodigal as input for diamond
    diamond_cmd += f"-d {db} "
    diamond_cmd += "-e 0.0001 " # maximum e-value of 10^-4
    diamond_cmd += "--id 35 " # minimum identity% coverage greater than or equal to 35%
    diamond_cmd += "--query-cover 50 " # query coverage greater than or equal to 50%
    diamond_cmd += "--subject-cover 50 " # reference coverage greater than or equal to 50%
    diamond_cmd += "-k 0 " # max number of target seq for align (default = 25)
    diamond_cmd += f"-p {threads} "
    diamond_cmd += f"-o {out}"
    
    # running the command
    p = subprocess.call(diamond_cmd, shell = True) # for testing printing the command
    
def parse_diamond(path):
    with open(path) as f:
        names = [
            "qseqid",
            "tseqid",
            "pid",
            "aln",
            "mis",
            "gap",
            "qstart",
            "qstop",
            "tstart",
            "tstop",
            "eval",
            "score",
        ]
        formats = [str, str, float, int, int, int, int, int, int, int, float, float]
        for line in f:
            values = line.split()
            yield dict([(names[i], formats[i](values[i])) for i in range(12)])

########################## Calculating PPQ/gPPQ Score ########################## 

def main():
    
    startTime = time.time()
    #run_prodigal(fna=args.fna, out=args.out)
    #run_diamond(faa=f'{args.out}/proteins.faa', db=f'{args.db}/proteins.dmnd', out=f'{args.out}/diamond.tsv', threads=args.threads)
    
    ########################## genome object ##########################

    ## collect genome information
    genomes = {}
    for header in SeqIO.parse(args.fna, 'fasta'):
        genome = Genome()
        genome.id = header.id.split()[0]
        genome.len = len(str(header.seq).upper())
        genomes[genome.id] = genome

    ########################## genes object ##########################
    
    genes = {}
    for header in SeqIO.parse(f'{args.out}/proteins.faa', 'fasta'):
        gene = Gene()
        gene.id = header.id.split()[0]
        gene.genome_id = header.id.rsplit('_',1)[0]
        gene.len = len(str(header.seq).upper())
        gene.ppq = None
        genes[gene.id] = gene
        genomes[gene.genome_id].genes.append(gene.id)


    ########################## reference ##########################
    
    reference = {}
    total_dbtype = []
    for row in csv.DictReader(open(f'{args.db}/proteins.tsv'), delimiter='\t'):
        ref = Reference()
        ref.gene_id = row['gene_id']
        ref.database = row['type']
        ref.genome_id = row['gene_id'].rsplit('_', 1)[0]
        reference[ref.gene_id] = ref
        total_dbtype.append(row['type'])

    totalvirus = total_dbtype.count('virus')
    totalplasmid = total_dbtype.count('plasmid')

    print(f'Total virus database count: {totalvirus} \nTotal plasmid database count: {totalplasmid}')


    ########################## PPQ/gPPQ ##########################
    
    ##### Store gene hits
    for row in parse_diamond(f'{args.out}/diamond.tsv'):
        # skip self hits
        if row['qseqid'].rsplit('_',1)[0] == row['tseqid'].rsplit('_',1)[0]:
            continue
        # plasmid hit
        elif reference[row['tseqid']].database == 'plasmid':
            genes[row['qseqid']].plasmid_hits.add(row['tseqid'].rsplit('_',1)[0])
        # virus hit
        elif reference[row['tseqid']].database == 'virus':
            genes[row['qseqid']].virus_hits.add(row['tseqid'].rsplit('_',1)[0])
    
    ##### Calculate PPQ scores
    for gene_id in genes:
        genes[gene_id].calc_ppq(totalvirus, totalplasmid)
        
    ##### Calculate gPPQ scores
    for genome_id in genomes:
        genomes[genome_id].calc_gppq(genes)

    ### writing gPPQ scores to a csv file ###
    with open(f'{args.out}/results.tsv', 'w') as f:
        fields = ['genome_id', 'num_genes', 'num_ppq', 'gppq', 'ppq_scores']
        f.write('\t'.join(fields)+'\n')
        for genome_id in genomes:
            row = {}
            row['genome_id'] = genome_id
            row['num_genes'] = genomes[genome_id].num_genes
            row['num_ppq'] = genomes[genome_id].num_ppq
            row['gppq'] = genomes[genome_id].gppq
            row['ppq_scores'] = ",".join(str(_) for _ in genomes[genome_id].ppq_scores)
            f.write('\t'.join(str(row[field]) for field in fields)+'\n')
        
if __name__ == "__main__":
    main()


