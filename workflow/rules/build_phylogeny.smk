import os
PATH = os.path.abspath('.')
RAXML_OUTPUT_FOLDER = f"{ PATH }/output/{ config['raxml_output_name'] }_RAxML"

include: 'utils.smk'
input_species = get_input_species()

rule fetch_sequences:
    output:
        sequences = 'output/odb_sequences.faa'
    params:
        all_species = input_species
    conda:
        '../envs/basic.yaml'
    log:
        'log/fetch_OrthoDB_sequences.log'
    script:
        '../scripts/fetch_odb_sequences.py'

'''
CAREFUL: if you want to provide a sequence file at this point,
make sure that each sequence is contained in one line, i.e. no
new line characters allowed inside sequences.
'''

checkpoint create_orthogroup_files:
    input:
        sequences = rules.fetch_sequences.output.sequences
    output:
        orthoGroups = directory('output/orthogroups/')
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/create_orthogroup_files.py'

# rule align:
#     input:
#         'output/orthogroups/{orthogroup}.fas'
#     output:
#         'output/alignments/{orthogroup}.aln'
#     conda:
#         '../envs/tree_building.yaml'
#     log:
#         'log/alignments/align_{orthogroup}.log'
#     shell:
#         'muscle -in {input} -out {output} -seqtype protein -quiet'

rule align:
    input:
        'output/orthogroups/{orthogroup}.faa'
    output:
        'output/alignments/{orthogroup}.aln'
    conda:
        '../envs/tree_building.yaml'
    log:
        'log/alignments/align_{orthogroup}.log'
    shell:
        'mafft --auto --leavegappyregion --anysymbol {input} > {output}'

rule trim:
    input:
        'output/alignments/{orthogroup}.aln'
    output:
        'output/alignments/{orthogroup}.aln.trm'
    conda:
        '../envs/tree_building.yaml'
    log:
        'log/trims/trim_{orthogroup}.log'
    shell:
        'trimal -in {input} -out {output} -strictplus'

def get_orthogroups(wildcards):
    orthogroups_dir = checkpoints.create_orthogroup_files.get(**wildcards).output[0] # get() here forces the checkpoint to rerun the DAG. E.g. without get(), I would only get a string of the output name.
    input = expand(rules.trim.output,orthogroup=glob_wildcards(os.path.join(orthogroups_dir, '{orthogroup}.fas')).orthogroup)
    return input

rule concatenate:
    input:
        get_orthogroups
    output:
        msa = 'output/multiple_sequence_alignment.faa'
    params:
        all_species = input_species
    conda:
        '../envs/tree_building.yaml'
    script:
        '../scripts/concatenate.py'

rule fasttree:
    input:
        rules.concatenate.output.msa
    output:
        'output/FastTree.tree'
    log:
        'log/FastTree.log'
    shell:
        'FastTree -lg -pseudo -gamma -spr 4 -mlacc 2 -slownni < {input} > {output}'

rule raxml:
    input:
        rules.concatenate.output.msa
    output:
        bestTree = f"{ RAXML_OUTPUT_FOLDER }/RAxML_bestTree.{ config['raxml_output_name'] }",
        bipartitions = f"{ RAXML_OUTPUT_FOLDER }/RAxML_bipartitions.{ config['raxml_output_name'] }",
        BranchLabels = f"{ RAXML_OUTPUT_FOLDER }/RAxML_bipartitionsBranchLabels.{ config['raxml_output_name'] }",
        bootstrap = f"{ RAXML_OUTPUT_FOLDER }/RAxML_bootstrap.{ config['raxml_output_name'] }",
        info = f"{ RAXML_OUTPUT_FOLDER }/RAxML_info.{ config['raxml_output_name'] }"
    conda:
        '../envs/tree_building.yaml'
    log:
        'log/RAxML.log'
    shell:
        "raxmlHPC-PTHREADS -T {config[raxml_cores]} -s {input} -n {config[raxml_output_name]} \
        -o {config[outgroup]} -N {config[bootstraps]} -m {config[model]} -f {config[algorithm]} \
        -x {config[rapid_bootstrap_seed]} -p {config[parsimony_seed]} -w {RAXML_OUTPUT_FOLDER}"
