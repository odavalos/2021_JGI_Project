#!/usr/bin/env python

import argparse
import pandas as pd

#### creating a text file to keep track of plasmid and phage proteins from the databases

# create a parser
parser = argparse.ArgumentParser()

# parser arguments
parser.add_argument('--input', type = argparse.FileType('r'), required = True, help = 'Input protein sequence file')
parser.add_argument('--db_type', type = str, required = True, help = 'Database type name to give new column to label header id\'s')
parser.add_argument('--merged_metadata', type = bool, default = False, help = 'whether extract metadata from both databases and to merge the files, default=False')
parser.add_argument('--input_2', type = argparse.FileType('r'), required = False, help = '2nd input sequence to merge the database metadata')
parser.add_argument('--db_type2', type = str, required = False, help = 'Second database type name to give new column to label header id\'s')
parser.add_argument('--output', type = str, required = True, help = 'Output file name')
parser.add_argument('--output_2', type = str, required = False, help = 'Second output file name')
parser.add_argument('--output_merged', type = str, required = False, help = 'Merged output file name')

args = parser.parse_args()


if args.merged_metadata == True:

	########################## Input_1 ##########################
	# read the fasta file and isolate header lines and separate the strings
	header_lines = []

	with args.input  as file:
		for line in file: 
			if line.startswith('>'):
				header_lines.append(line.strip())

			else:
				pass


	# clean header to only contain sequence names
	clean_header = []
	for line in header_lines:
		name = line.split(' ')[0]
		name = name[1:]
		clean_header.append(name)

	#print(f"First 10 names: {clean_header[:10]}")

	#print("#########################")
	#print("#########################")

	# clean names even further separate
	# Genbank Virus database input example: >[GENOME_ID]_[GENE_NUMBER] (ex: >AFYC01000010.1_3)

	genome_id = []
	gene_number = []

	for name in clean_header:
		if name.count('_') >= 2:
			mult_split = name.rsplit('_',1)
			genome_id.append(mult_split[0])
			gene_number.append(mult_split[1])

		else:
			splitname = name.split('_')
			genome_id.append(splitname[0])
			gene_number.append(splitname[1])

	#print(f"First ten genome id's: {genome_id[:10]} \nFirst ten gene_numbers:{gene_number[:10]}")

	db = [args.db_type] * len(genome_id)

	#create a dataframe that contains all information
	database_metainfo = pd.DataFrame({'sequence_header':clean_header, 'genome_id':genome_id, 'gene_number':gene_number, 'database':db})

	# save dataframe as csv file
	database_metainfo.to_csv(f'{args.output}_metadata.csv',index = False)
	
	print(f"Input_1 ==> First 10 names: {clean_header[:10]}")
	print("#########################")
	print("#########################")
	print(f"Input_1 ==> First ten genome id's: {genome_id[:10]} \nInput_1 ==> First ten gene_numbers:{gene_number[:10]}")
	print("#########################")
	print("#########################")
	print(f'Input_1 ==>First 10 rows of the csv:\n{database_metainfo[:10]}')

	########################## Input_2 ##########################
	# read the fasta file and isolate header lines and separate the strings
	header_lines = []
	with args.input_2  as file:
		for line in file: 
			if line.startswith('>'):
				header_lines.append(line.strip())

			else:
				pass


	# clean header to only contain sequence names
	clean_header = []
	for line in header_lines:
		name = line.split(' ')[0]
		name = name[1:]
		clean_header.append(name)

	#print(f"First 10 names: {clean_header[:10]}")

	#print("#########################")
	#print("#########################")

	# clean names even further separate
	# Genbank Virus database input example: >[GENOME_ID]_[GENE_NUMBER] (ex: >AFYC01000010.1_3)

	genome_id = []
	gene_number = []

	for name in clean_header:
		if name.count('_') >= 2:
			mult_split = name.rsplit('_',1)
			genome_id.append(mult_split[0])
			gene_number.append(mult_split[1])

		else:
			splitname = name.split('_')
			genome_id.append(splitname[0])
			gene_number.append(splitname[1])

	db = [args.db_type2] * len(genome_id)

	#create a dataframe that contains all information
	database_metainfo = pd.DataFrame({'sequence_header':clean_header, 'genome_id':genome_id, 'gene_number':gene_number, 'database':db})

	# save dataframe as csv file
	database_metainfo.to_csv(f'{args.output_2}_metadata.csv',index = False)

	print(f"Input_2 ==> First 10 names: {clean_header[:10]}")
	print("#########################")
	print("#########################")
	print(f"Input_2 ==> First ten genome id's: {genome_id[:10]} \nInput_2 ==> First ten gene_numbers:{gene_number[:10]}")
	print("#########################")
	print("#########################")
	print(f'Input_2 ==> First 10 rows of the csv:\n{database_metainfo[:10]}')

	########################## Merging Metadata ##########################

	db_1 = pd.read_csv(f'{args.output}_metadata.csv')
	db_2 = pd.read_csv(f'{args.output_2}_metadata.csv')
	merged_metadata = db_1.append(db_2, ignore_index=True) # does not attach second header as a row
	merged_metadata.to_csv('mergedDB_metadata.csv',index = False)


else:
	# read the fasta file and isolate header lines and separate the strings
	header_lines = []

	with args.input  as file:
		for line in file: 
			if line.startswith('>'):
				header_lines.append(line.strip())

			else:
				pass


	# clean header to only contain sequence names
	clean_header = []
	for line in header_lines:
		name = line.split(' ')[0]
		name = name[1:]
		clean_header.append(name)

	#print(f"First 10 names: {clean_header[:10]}")

	#print("#########################")
	#print("#########################")

	# clean names even further separate
	# Genbank Virus database input example: >[GENOME_ID]_[GENE_NUMBER] (ex: >AFYC01000010.1_3)

	genome_id = []
	gene_number = []

	for name in clean_header:
		if name.count('_') >= 2:
			mult_split = name.rsplit('_',1)
			genome_id.append(mult_split[0])
			gene_number.append(mult_split[1])

		else:
			splitname = name.split('_')
			genome_id.append(splitname[0])
			gene_number.append(splitname[1])

	#print(f"First ten genome id's: {genome_id[:10]} \nFirst ten gene_numbers:{gene_number[:10]}")

	db = [args.db_type] * len(genome_id)

	#create a dataframe that contains all information
	database_metainfo = pd.DataFrame({'sequence_header':clean_header, 'genome_id':genome_id, 'gene_number':gene_number, 'database':db})

	# save dataframe as csv file
	database_metainfo.to_csv(f'{args.output}_metadata.csv',index = False)

	print(f"Input_1 ==> First 10 names: {clean_header[:10]}")
	print("#########################")
	print("#########################")
	print(f"Input_1 ==> First ten genome id's: {genome_id[:10]} \nInput_1 ==> First ten gene_numbers:{gene_number[:10]}")
	print("#########################")
	print("#########################")
	print(f'Input_1 ==>First 10 rows of the csv:\n{database_metainfo[:10]}')