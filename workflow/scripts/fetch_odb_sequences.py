#!/usr/bin/env python3
## Loading the required packages
import sys
from collections import defaultdict
from datetime import datetime
import json
import requests

print('START: '+str(datetime.now())) ## Get starting time of script

host = snakemake.config["host"] ## Orthodb website that we will be using for querying

taxonomicLevel = snakemake.config["taxonomic_level"] ## taxonomic identifier (NCBI), e.g. ID for node Hymenoptera

## Dictionary of all species, pairing numerical identifiers to a human readable identifier
speciesDict = snakemake.params.species

numberOfSpecies = len(speciesDict) ## The number of species we are looking at

limit = snakemake.config["request_limit"] ## A limit for getting results from the queries we will be submitting
universalPercentage = snakemake.config["universality"] ## The percentage of species required to have this gene in order for the query to return a result
singleCopyPercentage = snakemake.config["singlecopyness"] ## The percentage of number of species that have this gene as a single opy

print('\nuniversalPercentage: ' + universalPercentage)
print('singleCopyPercentage: ' + singleCopyPercentage + '\n')

query = 'search?limit='+limit+'&level='+taxonomicLevel+'&universal='+universalPercentage+'&singlecopy='+singleCopyPercentage ## Our query
print('Your query: ' +host+query+'\n')
response = requests.get(host+query) ## Gets a response from the querying the url (host + query)
print('Response code: '+str(response.status_code)+'\n')

while response.status_code == 504: ## 504 Gateway Timeout error is an HTTP status code that means that one server did not receive a timely response from another server. So we try again until the query works.
    response = requests.get(host+query)
    print('Gateway Timeout error. Response code: '+str(response.status_code)+'\n')
    print('Trying again ...')
    # print('Retrying in 5 seconds ...')
    # time.sleep(5)
if response.status_code != 200: ## Any other response code than 200 is an error
    print('Bad query: '+host+query+'\n')
    print('Response code: '+str(response.status_code)+'\n')
    print('Response content: '+response.text+'\n')
    sys.exit() ## Quits the script, because we have a bad query

jsonData = json.loads(response.content.decode('utf-8')) ## Convert the response we got from the api to python json format
orthoGroups = jsonData['data'] ## Get all the orthologous groups from the JSON variable

print('Your query returned '+str(len(orthoGroups))+' orthologous groups'+'\n')
print('Now retrieving sequences for '+str(numberOfSpecies)+' species: '+', '.join(speciesDict.keys())+'\n')

outfile = open(snakemake.output["sequences"],'w') ## We will write the sequences that are single copy and universal to this fasta file
speciesList = ','.join(speciesDict.keys()) ## Take all keys (= numerical species identifiers) from our species dictionary and convert it to a comma separated list
processedGroupCount = 0
usedGroupCount = 0
present_in_all_count = 0
single_copy_count = 0

for orthoGroup in orthoGroups: ## Query over each orthologous group
    speciesCopyDict = defaultdict(list) ## A dictionary where we will store all genes of the orthologous group belonging to a species
    query = 'fasta?id='+orthoGroup+'&species='+speciesList ## Our query
    response = requests.get(host+query) ## Getting the response from the API url (host + query)
    if response.status_code != 200: ## Any response code thats not 200 means we have a bad query
        print('Bad query: '+host+query+'\n')
        print('Response code: '+str(response.status_code)+'\n')
        print('Response content: '+response.text+'\n')
        sys.exit() ## Exit the script, because we have a bad query
    else:
        fasta = response.text ## Each response consists of fasta files
        for line in fasta.split('\n'): ## Split the fasta variabe on lines (= \n) to s we can iterate over each line
            if line.startswith('>'): ## If the line starts with '>' it is the header of the fasta file
                speciesIdentifier = line.strip('>').split(':')[0] ## We isolate the numerical species identifier
                speciesName = speciesDict[speciesIdentifier] ## We convert the species numerical identifier to the name we gave it in our species dictionary at the start of this script
                geneIdentifier = line.split('"pub_gene_id":"')[1].split('"')[0] ## We isolate the gene identifier
            elif not line.strip() == '': ## Anything that doesn't start with '>' is the protein sequence
                sequence = line.strip() ## To remove the line break at the end of the line
                speciesCopyDict[speciesName].append(speciesName+'|'+orthoGroup+'|'+geneIdentifier+'\n'+sequence) ## Each fasta file is stored in the species copy dict, for each species, in list format
    isUniversal = True ## A boolean value to check if we want to write the sequences for this group to our outfile
    isSingleCopy = True

    ## Code block that selects only single-copy orthogroups and present-in-all species. Against the purpose of user selection of query paramters.
    if not len(speciesCopyDict.keys()) == numberOfSpecies: ## If the amount of species for which we have an orthologous gene doesnt match all the species we are looking for, change the useGroup boolean to False because the orthologous group is not universal
        isUniversal = False
    # else: ## If the orthologous group is not single copy in every species, change the useGroup to false
    for species,genes in speciesCopyDict.items():
        if not len(genes) == 1:
            isSingleCopy = False

    if isUniversal == True and isSingleCopy == True: ## If the useGroup is still True, write the sequences to our outfile
        for species,genes in speciesCopyDict.items():
            for gene in genes:
                outfile.write('>'+gene+'\n')
        usedGroupCount += 1
        present_in_all_count += 1
        single_copy_count += 1
    elif isUniversal == True and isSingleCopy == False:
        for species,genes in speciesCopyDict.items():
            gene_lengths = {}
            for gene in genes:
                gene_id = gene.split('|')[2]
                gene_seq = gene.split('\n')[1]
                gene_length = len(gene_seq)
                gene_lengths[gene_length] = gene
            max_length = max(gene_lengths.keys())
            outfile.write('>'+gene_lengths[max_length]+'\n')
        usedGroupCount += 1
        present_in_all_count += 1
    elif isUniversal == False and isSingleCopy == True:
        for species,genes in speciesCopyDict.items():
            for gene in genes:
                outfile.write('>'+gene+'\n')
        last_gene = gene.split('\n')[1]
        gene_length = len(last_gene)
        for species in speciesDict.keys():
            if species not in speciesCopyDict.keys():
                outfile.write('>'+species+'|'+orthoGroup+'|'+'gap'+'\n'+'-'*gene_length+'\n')
        usedGroupCount += 1
        single_copy_count += 1
    elif isUniversal == False and isSingleCopy == False:
        for species,genes in speciesCopyDict.items():
            gene_lengths = {}
            for gene in genes:
                gene_id = gene.split('|')[2]
                gene_seq = gene.split('\n')[1]
                gene_length = len(gene_seq)
                gene_lengths[gene_length] = gene
            max_length = max(gene_lengths.keys())
            outfile.write('>'+gene_lengths[max_length]+'\n')
        for species in speciesDict.keys():
            if species not in speciesCopyDict.keys():
                outfile.write('>'+species+'|'+orthoGroup+'|'+'gap'+'\n'+'-'*max_length+'\n')
        usedGroupCount += 1

    processedGroupCount += 1
    if processedGroupCount % 10 == 0: ## Progress printout to terminal
        print(str(processedGroupCount)+' orthologous groups processed... '+str(datetime.now()))

outfile.close()
print('Processed '+str(usedGroupCount)+' orthologous groups across '+str(numberOfSpecies)+' species\n')
print('Of which: ' + str(single_copy_count) ' single-copy and ' + str(present_in_all_count) + 'present in all species.')
print('END: '+str(datetime.now()))
