import os

absolute = os.path.abspath(os.getcwd())
absolute = absolute + '/../'

data = {
    "visualizing": absolute + 'data/visualizing',
    "rearrange": absolute + 'data/rearrange',
    "sequences_regions": absolute + 'data/sequences_regions',
    "clustering": absolute + 'data/clustering',
    "filtering": absolute + 'data/filter',
    "alignments": absolute + 'data/alignments_global_needle',
}

fasta = {
    "sequences": absolute + 'data/disordered_sequences_splitted_fasta',
    "regions": absolute + 'data/disordered_regions_splitted_fasta',
}

alignments = {
    "global": absolute + 'data/alignments_global_needle',
}


bash = {
    "results": absolute + 'bash_scripts/results_needle',
}
