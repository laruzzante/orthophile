# orthogeny

Snakemake pipeline to build a tree from OrthoDB single-copy orthologues. 


# usage:

Install conda or miniconda (https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

Download the orthogeny git repository.

Open a terminal and move into the <orthogeny> folder:
`cd <your/path/to/orthogeny>`

Create the specific conda environment:
`conda env create -f environment.yaml`

Activate the conda environment:
`conda activate orthogeny`

Move into the <workflow> folder:
`cd workflow`

Edit the configuration file (specify yout OrthoDB species list, orthologous group selection parameters, RAxML parameters):
config.yaml

Rune the pipeline, where <N> is the number of cores you want snakemake to use:
`snakemake --cores <N>`

At the end of the computation, you will see a Snakemake error message but the RAxML results will still appear at the current location.
