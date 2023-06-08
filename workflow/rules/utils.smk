### This function reads the input species list and creates a species dictionary with the species taxonomy IDs (taxid) as keys,
### and the species name as values. Careful about the _0 present in OrthoDB nomenclature but absent in NCBI Taxonomy ids.

def get_input_species():
    species = {}
    with open(config['input_species']) as f:
        for line in f:
            if(',') in line:
                name = line.strip().split(', ')[1]
                if('_0') not in line: ## _X in OrthoDB species IDs indicates that a species may have multiple indexed genomes. In Arthropods, there is usually only 1 genome per species
                ## thus we usually only find _0. This may change with species in other lineages, f.eg. yeasts and bacteria. NCBI taxids do not have an underscore after the value.
                    id = line.strip().split(', ')[0] + '_0' ## If species ID do not have the _0 ending typical of OrthoDB arthropods, then we add it
                else:
                    id = line.strip().split(', ')[0]
            else:
                id = line.strip().split(', ')[0]
                name = id ## If no species name is provided in the input species list, then we set it as the name as their provided NCBI tax id
            species[id] = name
    return species
