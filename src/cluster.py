#!/usr/bin/env python3

import argparse
import os
import sys
import networkx as nx
import pandas as pd
import config as cfg


def main():
    cluster_file()
    # filter_cluster()

# remove uniprots that do not have a structural state region associated

def filter_cluster():

    removedproteins = 0
    regionsarr = []
    uniprots = []
    with open(cfg.data['sequences_regions'] + '/search-disprot-filtered.tsv', 'r') as regions:

        for line in regions:
            region = line.split('\t')[5]
            uniprot = line.split('\t')[0]
            uniprots.append(uniprot)
            regionsarr.append(region)
    # print(regionsarr)
    with open(cfg.data['clustering'] + '/components-filtered.tsv', 'w+') as outcomponents:
        outcomponents.write('reg1' + '\t' + 'reg2' + '\t' + 'id1' + '\t' + 'id2' + '\t'
                            + 'l1' + '\t' + 'l2' + '\t' + 'identity' + '\t' + 'similarity' + '\t'
                            + 'gaps' + '\t' + 'score' + '\t'
                            + 'align1' + '\t' + 'align2' + '\t'

                            + 'regions1' + '\t' + 'regions2' + '\t'
                            + 'cluster' + '\n')
        with open(cfg.data['clustering'] + '/components.tsv', 'r') as alignments:
            for alignment in alignments:

                id1 = alignment.split('\t')[2]
                id2 = alignment.split('\t')[3]
                reg1 = alignment.split('\t')[11]
                reg2 = alignment.split('\t')[12]
                filteredreg1 = []
                filteredreg2 = []

                if id1 in uniprots and id2 in uniprots:

                    for reg in reg1.split(','):
                        if reg.split('_')[0] in regionsarr:
                            filteredreg1.append(reg)
                    for reg in reg2.split(','):
                        if reg.split('_')[0] in regionsarr:
                            filteredreg2.append(reg)

                    outcomponents.write(alignment.split('\t')[0] + '\t' + alignment.split('\t')[1] + '\t' +
                                        alignment.split('\t')[2] + '\t' + alignment.split('\t')[3] + '\t'
                                        + alignment.split('\t')[4] + '\t' + alignment.split('\t')[5] + '\t'
                                        + alignment.split('\t')[6] + '\t' + alignment.split('\t')[7] + '\t'
                                        + alignment.split('\t')[8] + '\t' + alignment.split('\t')[9] + '\t'
                                        + alignment.split('\t')[10] + '\t' + '\t'

                                        + ','.join(filteredreg1) + '\t'
                                        + ','.join(filteredreg2) + '\t'
                                        + alignment.split('\t')[13])
                else:
                    removedproteins += 1
    print('Removed after filtering: ' + str(removedproteins))


# to be executed after having filtered the results
# generate clusters
def cluster_file():

    index = []
    idOne = []
    idTwo = []
    with open(cfg.data['filtering'] + '/filter-needle-identity.txt', 'r') as input:
        count = 0

        for row in input:
            idOne.append(row.split('\t')[2])
            idTwo.append(row.split('\t')[3])
            index.append(count)
            count += 1

    d = {"id1": pd.Series(idOne, index=index), "id2": pd.Series(idTwo, index=index)}
    df = pd.DataFrame(d)
    G = nx.from_pandas_edgelist(df, 'id1', 'id2')

    components = {}
    count = 0
    for comp in list(nx.connected_components(G)):
        components[count] = comp
        count += 1

    componentsDict = {}
    for key in components:
        for value in components[key]:
            componentsDict[value] = key
    print(componentsDict)

    with open(cfg.data['filtering'] + '/filter-needle-identity.txt', 'r') as table:
        with open(cfg.data['clustering'] + '/components.tsv', 'w+') as tablecomp:
            for line in table:
                id = line.split('\t')[2]
                compVal = componentsDict[id]
                print(compVal)
                tablecomp.write(line.strip() + '\t' + str(compVal) + '\n')


def parse_arguments(initargs=None):
    if initargs is None:
        initargs = []
    parser = argparse.ArgumentParser(description="Parse input alignment file.")
    parser.add_argument('-i', help="Input alignment file")
    if not initargs or len(sys.argv) > 1:
        args = parser.parse_args()
    else:
        args = parser.parse_args(initargs)
    return args


if __name__ == '__main__':
    sys.exit(main())
