#!/usr/bin/env python3
import sys
from datetime import datetime
import os
import shutil

print('START: '+str(datetime.now())+'\n') ## Starting time of the analysis

infile = open(snakemake.input["sequences"])
outdir = snakemake.output[0]

if not os.path.isdir(outdir):
    os.mkdir(outdir)

for line in infile.readlines(): ## Iterate over the input file
	if line.startswith('>'): ## Take the species and orthology group from the header line, create new header line with just the species name
		ssline = line.strip().split('|')
		species = ssline[0]
		orthoGroup = ssline[1]
		header = species
	else: ## If line doesn't start with '>' its the sequence
		sequence = line
		outfile = open(os.path.join(outdir,orthoGroup+'.fas'),'a') ## Write the fasta to the orthoGroup file
		outfile.write(header+'\n'+sequence) ## Header consists of species name only
		outfile.close()

infile.close()
