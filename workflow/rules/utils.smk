def get_input_species():
    species = {}
    with open(config['input_species']) as f:
        for line in f:
            id = line.strip().split(',')[1] + '_0'
            name = line.strip().split(',')[0]
            species[id] = name
        return species
