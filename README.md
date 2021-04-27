# Orthophile

Snakemake pipeline to build phylogenetic trees with thousand of orthologues genes from ***[OrthoDB](https://www.orthodb.org/)***.


## Usage

Install conda or miniconda (https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

Download the orthophile git repository.

Open a terminal and move into the **orthophile** folder:
`cd <your/path/to/orthophile>`

Move into the **workflow** folder:
`cd workflow`

Prepare an input file containing the list of species that you would like to include in the phylogeny, using this format (one species per row):
`species_name,NCBI_taxonomy_ID`
where 'species_name' can be any name you want to associate to that species, followed by a comma and its corresponding NCBI Taxonomy identifier.

Edit the configuration file **config.yaml** (specify your orthologous group selection parameters and RAxML parameters).

Rune the pipeline, where **N** is the number of cores you want snakemake to use:
`snakemake --cores <N> --use-conda`

The RAxML results will be stored in a specific folder in the output.

