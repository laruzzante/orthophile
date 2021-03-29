#!/usr/bin/env python3
from collections import defaultdict

input = snakemake.input
all_species = snakemake.params.all_species

## We're concatenating the alignments into one big alignment, as this is more useful for a phylogenetic tree than thousands of different alignments

concatenatedAlignmentDictionary = defaultdict(str) ## Dictionary to store the concatenated alignments
for fin in input:
    orthogroup_species_list = []
    for line in open(fin):
        if line.startswith('>'):
            seq_length = 0
            species = line[1:].strip()
            orthogroup_species_list.append(species)
        else:
            sequence_part = line.strip()
            concatenatedAlignmentDictionary[species] += sequence_part
            seq_length += len(sequence_part)
    for species_id,species in all_species.items():
        if species not in orthogroup_species_list:
            print(f'\n{speceis} missing from the alignment file {fin}, adding gap sequence of same length.')
            if species not in concatenatedAlignmentDictionary.keys():
                concatenatedAlignmentDictionary[species] = '-'*seq_length+'\n'
            else:
                concatenatedAlignmentDictionary[species] += '-'*seq_length

outfile = open(snakemake.output["msa"],'w')
for species,sequence in concatenatedAlignmentDictionary.items(): ## Every species with its concatenated sequence is written out to the concatenated alignment file
	outfile.write('>'+species+'\n')
	outfile.write(sequence+'\n')
outfile.close()
