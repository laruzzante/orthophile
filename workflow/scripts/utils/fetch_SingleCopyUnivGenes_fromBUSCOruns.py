#!/usr/bin/env python3

'''
    Author: Livio Ruzzante

    Script that extracts the BUSCO single copy ohrtologs present across
    all species.
'''

import os
from collections import defaultdict


path = '/home/lruzzant/Phylogenies/Genomes' + \
       '/busco_results/MosquitoPhylogenies-project/'
busco_runs = os.listdir(path)

# Hardcoded deletion of An. maculatus BUSCO run (bad assembly) to obtain 2x more genes.
try:
    deleted_run = "An_maculatus"
    busco_runs.remove(deleted_run)
except ValueError:
    print("ERROR: run not found, could not delete", deleted_run)

buscoHitsDict = defaultdict(list)
for run in busco_runs:
    sc_path = os.path.join(path, run, 'single_copy_busco_sequences/')
    busco_hits = os.listdir(sc_path)
    # Keeping only .faa files (protein sequence and not nucleotide sequence)
    busco_hits_faa = []
    for i in range(0, len(busco_hits)):
        if busco_hits[i].endswith('.faa'):
            busco_hits_faa.append(busco_hits[i])
    buscoHitsDict[run] = busco_hits_faa


# Take the intersection of lists of orthologs to end up with
# a set of common orthologs.

# Manually selecting a first BUSCO run to start with:
first_run = 'Ae_aegypti'
buscoHits_intersection = buscoHitsDict[first_run]
for run in buscoHitsDict.keys():
    buscoHits_intersection = list(set(buscoHits_intersection).
                                  intersection(buscoHitsDict[run]))

print('\nUniversal single copy genes found with BUSCO across',
      len(buscoHitsDict.keys()), 'runs:', len(buscoHits_intersection))
print('\nIncluded runs:')
print('\n'.join(sorted(busco_runs)))

# Collating all UN SC gene sequences from all species into a single file
allSeqDict = defaultdict(lambda: defaultdict())
for og in sorted(buscoHits_intersection):
    for run in sorted(buscoHitsDict.keys()):
        # Checking if every species has all the SC orthogroups
        # of the intersection.
        if og not in buscoHitsDict[run]:
            print('ERROR: ', og, 'missing in ', run)
        else:
            with open(path + run + '/single_copy_busco_sequences/' +
                      og, 'r') as f:
                lines = f.readlines()
                og_id = og.split('.')[0]
                allSeqDict[og_id][run] = lines[1].strip()

# Printing to file
with open('BUSCO_hits_sequences.faa', 'w') as f:
    for og_id in sorted(allSeqDict.keys()):
        for run in sorted(allSeqDict[og_id].keys()):
            f.write('>' + run + '|' + og_id + '\n')
            f.write(allSeqDict[og_id][run] + '\n')
