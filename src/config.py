import os

absolute = os.path.abspath(os.getcwd())
absolute = absolute + '/../'

data = {
    "visualizing": absolute + 'data/visualizing',
    "rearrange": absolute + 'data/rearrange',
    "sequences_regions": absolute + 'data/sequences_regions',
    "clustering": absolute + 'data/clustering',
    "filtering": absolute + 'data/filter',
}

fasta = {
    "sequences": absolute + 'data/disordered_sequences_splitted_fasta',
    "regions": absolute + 'data/disordered_regions_splitted_fasta',
}

alignments = {
    "global": absolute + 'data/alignments_global_needle',
}