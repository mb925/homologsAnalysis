import config as cfg
import pandas as pd
import sys
import math


def main():
    calculate_regions()
    # calculate_sequences()
    # calculate_regions_size()


def calculate_regions_size():
    regs = pd.read_csv(cfg.data['sequences_regions'] + '/search_in_disprot.tsv', sep='\t').drop_duplicates(
        subset=['disprot_id'], keep=False)
    print(regs)


def calculate_sequences():
    # added temporary header to calculate this
    ids1 = pd.read_csv(cfg.data['filtering'] + '/filter-needle-identity.txt', sep='\t')['id1']

    ids2 = pd.read_csv(cfg.data['filtering'] + '/filter-needle-identity.txt', sep='\t')['id2']

    ids = pd.concat([ids1, ids2]).drop_duplicates()
    df = pd.DataFrame()
    df['disprot_id'] = ids
    calculate_regions()


def calculate_regions():
    # ids1 = pd.read_csv(cfg.data['alignments'] + '/global_triangular_matrix.csv', sep='\t')['id1']
    #
    # ids2 = pd.read_csv(cfg.data['alignments'] + '/global_triangular_matrix.csv', sep='\t')['id2']
    #
    # ids = pd.concat([ids1, ids2]).drop_duplicates()
    # df = pd.DataFrame()
    # df['acc'] = ids

    # ids1 = pd.read_csv(cfg.data['filtering'] + '/filter-needle-identity.txt', sep='\t')['id1']
    #
    # ids2 = pd.read_csv(cfg.data['filtering'] + '/filter-needle-identity.txt', sep='\t')['id2']
    #
    # ids = pd.concat([ids1, ids2]).drop_duplicates()
    # df = pd.DataFrame()
    # df['disprot_id'] = ids

    # reg = pd.read_csv(cfg.data['sequences_regions'] + '/search_in_disprot.tsv', sep='\t')
    # merged_df = reg.merge(df, on='disprot_id')
    # merged_df = reg.drop_duplicates(subset=['acc', 'start', 'end'], keep=False)
    # size = len(merged_df)
    # tot_reg = pd.read_csv(cfg.bash['results'] + '/all-global-needle.txt', sep='\t')['reg2'].dropna()
    tot_reg = pd.read_csv(cfg.data['filtering'] + '/filter-needle-identity.txt', sep='\t')['reg1'].dropna()

    count = 0
    for row in tot_reg:
        print(row)

        values = row.split(',')
        for value in values:
            count += 1
            print(count)


if __name__ == '__main__':
    sys.exit(main())
