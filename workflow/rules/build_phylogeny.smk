import os
PATH = os.path.abspath('.')
RAXML_OUTPUT_FOLDER = f"{ PATH }/output/{ config['raxml_output_name'] }_RAxML"

rule fetch_sequences:
    output:
        sequences = 'output/odb_sequences.fasta'
    conda:
        '../envs/basic.yaml'
    script:
        'scripts/fetch_odb_sequences.py'

rule align_trim_concatenate:
    input:
        sequences = rules.fetch_sequences.output.sequences
    output:
        msa = 'output/multiple_sequence_alignment.fasta'
    conda:
        '../envs/tree_building.yaml'
    script:
        'scripts/align_trim_concatenate.py'

rule build_tree:
    input:
        rules.align_trim_concatenate.output.msa
    output:
        bestTree = f"{ RAXML_OUTPUT_FOLDER }/RAxML_bestTree.{ config['raxml_output_name'] }",
        bipartitions = f"{ RAXML_OUTPUT_FOLDER }/RAxML_bipartitions.{ config['raxml_output_name'] }",
        BranchLabels = f"{ RAXML_OUTPUT_FOLDER }/RAxML_bipartitionsBranchLabels.{ config['raxml_output_name'] }",
        bootstrap = f"{ RAXML_OUTPUT_FOLDER }/RAxML_bootstrap.{ config['raxml_output_name'] }",
        info = f"{ RAXML_OUTPUT_FOLDER }/RAxML_info.{ config['raxml_output_name'] }"
    conda:
        '../envs/tree_building.yaml'
    log:
        'log/RAxML_log.txt'
    shell:
        "raxmlHPC-PTHREADS -T {config[raxml_cores]} -s {input} -n {config[raxml_output_name]} \
        -o {config[outgroup]} -N {config[bootstraps]} -m {config[model]} -f {config[algorithm]} \
        -x {config[rapid_bootstrap_seed]} -p {config[parsimony_seed]} -w {RAXML_OUTPUT_FOLDER}"
