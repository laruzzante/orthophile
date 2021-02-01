#!/usr/bin/env python3
import sys
from datetime import datetime
from collections import defaultdict
import subprocess
import os
import shutil

print('START: '+str(datetime.now())+'\n') ## Starting time of the analysis

infile = open(snakemake.input["sequences"])

relative_path = ''

## Create a clean directory to store the orthology groups
if not os.path.exists(relative_path + 'orthoGroups/'):
	os.makedirs(relative_path + 'orthoGroups/')
else:
	shutil.rmtree(relative_path + 'orthoGroups/')
	os.makedirs(relative_path + 'orthoGroups/')

for line in infile.readlines(): ## Iterate over the input file
	if line.startswith('>'): ## Take the species and orthology group from the header line, create new header line with just the species name
		ssline = line.strip().split('|')
		species = ssline[0]
		orthoGroup = ssline[1]
		header = species
	else: ## If line doesn't start with '>' its the sequence
		sequence = line
		outfile = open(relative_path + 'orthoGroups/'+orthoGroup+'.fas','a') ## Write the fasta to the orthoGroup file
		outfile.write(header+'\n'+sequence) ## Header consists of species name only
		outfile.close()

infile.close()

## Create a clean alignments directory where we store the alignments
if not os.path.exists(relative_path + 'alignments/'):
	os.makedirs(relative_path + 'alignments/')
else:
	shutil.rmtree(relative_path + 'alignments/')
	os.makedirs(relative_path + 'alignments/')

for fn in os.listdir(relative_path + 'orthoGroups/'): ## Iterate over the orthology group fasta files, run the muscle alignment program
	if fn.endswith('.fas'):
		alnFile = fn.replace('.fas','.aln')
		subprocess.Popen('muscle -in '+ relative_path + 'orthoGroups/'+fn+' -out '+ relative_path + 'alignments/'+alnFile+' -quiet' ,shell=True)

for fn in os.listdir(relative_path + 'alignments/'): ## Iterate over the alignment files for each orthology group, run the trimal program to remove non-conserved regions to improve the alignment
	if fn.endswith('.aln'):
		alnFile = fn
		trmFile = fn.replace('.aln','.trm.aln')
		subprocess.Popen('trimal -in ' + relative_path + 'alignments/'+alnFile+' -out ' + relative_path + 'alignments/'+trmFile+' -strictplus', shell=True)

## We're concatenating the alignments into one big alignment, as this is more useful for a phylogenetic tree than 1500 different alignments
concatenatedAlignmentDictionary = defaultdict(str) ## Dictionary to store the concatenated alignments
for fn in os.listdir(relative_path + 'alignments/'): ## Iterating over all alignments
		if fn.endswith('.trm.aln'):
			species = ''
			for line in open(relative_path + 'alignments/'+fn): ## Every fasta sequence for every alignment file is appended to a dictionary with a 'species' -> 'sequence' pair
				if line.startswith('>'):
					if not species == '':
						concatenatedAlignmentDictionary[species] += sequence
					species = line.strip('>').strip().split(' ')[0]
					sequence = ''
				else:
					sequence += line.strip()
			concatenatedAlignmentDictionary[species] += sequence

outfile = open(snakemake.output["msa"],'w')
for species,sequence in concatenatedAlignmentDictionary.items(): ## Every species with its concatenated sequence is written out to the concatenated alignment file
	outfile.write('>'+species+'\n')
	outfile.write(sequence+'\n')
outfile.close()
