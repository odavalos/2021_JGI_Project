#!/usr/bin/env python

# std libraries
import os
import sys
import subprocess
import argparse
import pandas as pd

# creating the parser
parser = argparse.ArgumentParser()

# parser arguments
parser.add_argument('--fna', type = str, required = True, help = 'Input fasta file or path to the file')
parser.add_argument('--db', type = str, required = True, help = 'Database file or path to database')
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
    2. Output as genbank format (might change to fna?)
    """
    # get filename from path
    base = os.path.basename(f'{input}') # currently have set to zero since the directory contains two files
    base = os.path.splitext(base)[0] # the output is a tuple that where the first value the name
    
    # creating the shell script command
    prod_cmd = "prodigal "
    prod_cmd += f"-i {input} " # input is fna 
    prod_cmd += f"-o {base}_prodigal.gbk " # can it write out to faa?
    prod_cmd += f"-a {base}_prodigal.faa " # this is an output # give same name
    prod_cmd += "-p meta"
    
    # running the command
    p = subprocess.call(prod_cmd, shell = True) # for testing printing the command
    # add subprocess.stderr=PIPE
    # out, err = p.communicate()
    # print(err)


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

def run_diamond(input_faa, db, output):
    """
    Running diamond
    """
    
    # creating the shell script with the command
    diamond_cmd = "diamond "
    diamond_cmd += "blastp "
    diamond_cmd += f"-q {input_faa} " # use output from prodigal as input for diamond
    diamond_cmd += f"-d {db} "
    diamond_cmd += "-e 0.0000099 " # maximum e-value of 10^-4
    diamond_cmd += "--id 35 " # minimum identity% coverage greater than or equal to 50%
    diamond_cmd += "--query-cover 50 " # query coverage greater than or equal to 50%
    diamond_cmd += "--subject-cover 50 " # reference coverage greater than or equal to 50%
    diamond_cmd += f"-o {output}_dmnd.tsv"
    
    # running the command
    p = subprocess.call(diamond_cmd, shell = True) # for testing printing the command

########################## Calculating PPQ/gPPQ Score ########################## 

# create a dataframe with pandas
diamond_df = pd.read_csv(f'{inputbase}_dmnd.tsv', sep='\t', header = None)

# list of column names
output_colnames_dmnd = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

# change the dataframe column names
diamond_df.columns = output_colnames_dmnd

# read in merged meta data csv file
mergedDB_metadata = pd.read_csv('mergedDB_metadata.csv')

# match hits with meta data
hits = pd.merge(diamond_df,mergedDB_metadata, left_on = 'sseqid', right_on = 'sequence_header')


virus_totaldb = mergedDB_metadata.query('database == "virus"')['database'].count()
plasmid_totaldb = mergedDB_metadata.query('database == "plasmid"')['database'].count()

print(f"The total virus database is: {virus_totaldb}")
print(f"The total plasmid database is: {plasmid_totaldb}")

# count all the matches
matches = hits.groupby(['qseqid', 'database'])['database'].size().reset_index(name='count')


def calc_ppq(df):
    """
    Function calculates the qPPQ scores from Pfeifer et.al. 2021.
    1. creates empty lists for inputs
    2. parses the data
        - multi_matches: sequences that contain hits for both phages and plasmids
        - unique_matches: sequences that match up only to phages or plasmids
    3. loops through the parsed dataframes and calculate the ppq for each. 
        - calculates ppq for multi_matches first
        - calculates ppq for multi_matches second
    4. merges appended sequence_id's and ppq_scores into a dataframe as the output
    """
    sequence_id = []
    ppq_score = []
    
    
    # split the dataframe i.e. between duplicates (sequences that match both virus and plasmid) and non duplicates
    
    # virus and plasmid matching sequences
    multi_matches = df[df.duplicated('qseqid', keep = False)] # duplicates i.e. matches that contain virus and plasmid matches
    
    # uniquely matching sequences i.e. virus or plasmid
    unique_matches = test_matches.drop_duplicates(subset = ['qseqid'], keep = False) # all unique ids
    
    ####### calculating ppq for multi matching sequences #######
    mm_unique = multi_matches[['qseqid']].drop_duplicates()['qseqid'].tolist()
    for i in mm_unique:
        sequence_id.append(i)
        mm_unique_sub = multi_matches[multi_matches['qseqid'] == i]
        virus_num = int(mm_unique_sub[mm_unique_sub['database'] == 'virus']['count'].values)
        plasmid_num = int(mm_unique_sub[mm_unique_sub['database'] == 'plasmid']['count'].values)
        ppq = (virus_num / virus_totaldb) / ((virus_num/virus_totaldb)+(plasmid_num/plasmid_totaldb))
        ppq_score.append(ppq)
    
    ####### calculating ppq for unique matching sequences #######
    for index, row in unique_matches.iterrows():
        sequence_id.append(row['qseqid'])
        if row['database'] == 'virus':
            ppq = (row['count']/virus_totaldb) / ((row['count']/virus_totaldb)+(0/plasmid_totaldb))
            ppq_score.append(ppq)
        else:
            ppq = (0/virus_totaldb) / ((0/virus_totaldb)+(row['count']/plasmid_totaldb))
            ppq_score.append(ppq)
    
    # report results as dataframe
    ppq_df = pd.DataFrame({'Sequence_ID':sequence_id, 'PPQ_Score':ppq_score})

    # save dataframe as csv file
    ppq_df.to_csv(f'{inputbase}_ppqscores.csv',index = False)

    return(ppq_df)


########################## Running the program ########################## 
run_prodigal(args.fna)
make_diamonddb(args.db)
run_diamond(input_faa=f'{inputbase}_prodigal.faa', db=f'{db_base}.dmnd', output=f'{inputbase}')
calc_ppq(matches)



