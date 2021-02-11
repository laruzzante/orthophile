# Orthophile

Snakemake pipeline to build phylogenetic trees from ***[OrthoDB](https://www.orthodb.org/)*** universal single-copy orthologues.


## Usage

Install conda or miniconda (https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

Download the orthophile git repository.

Open a terminal and move into the **orthophile** folder:
`cd <your/path/to/orthogeny>`

Create the specific conda environment:
`conda env create -f environment.yaml`

Activate the conda environment:
`conda activate orthophile`

Move into the **workflow** folder:
`cd workflow`

Edit the configuration file **config.yaml** (specify yout OrthoDB species list, orthologous group selection parameters and RAxML parameters).

Rune the pipeline, where **N** is the number of cores you want snakemake to use:
`snakemake --cores <N>`

At the end of the computation, you will see a Snakemake error message but the RAxML results will still appear at the current location.

Once you are done with Orthophile, you can deactivate the conda environment by either closing the terminal or by executing:
`conda deactivate`
