#!/usr/bin/env python3
import config as cfg
import argparse
import os
import sys
import re

def main():
    parse_file()

#regions
#That info is not really important because downstream I will substitute those regions with the converted ones
# it is probably not useful to add the region info already
# I probably could filter here already for identity and structural state, but I will not know if needle calculates everything
def parse_file():
    args = parse_arguments()

    id1 = ''
    id2 = ''
    seq1 = ''
    seq2 = ''
    len1 = ''
    reg1 = []
    reg2 = []
    disprotid1 = ''
    disprotid2 = ''
    uniprot = ''
    region = ''
    start = ''
    end = ''
    count = 0
    flag = False
    alignments = ''
    if args.i:
        inf = open(args.i, 'r')
    else:
        inf = sys.stdin


    # regexpUniProt = re.compile('^[OPQopq][0-9][A-Za-z0-9]{3}[0-9]|[A-Na-nR-Zr-z][0-9]([A-Za-z][A-Za-z0-9]{2}[0-9]){1,2}$')  # seqres
    for line in inf:
        # print(line)

        if line.startswith("# 1:"):
            id1 = line.split(":")[1].strip()
            # if regexpUniProt.match(id1) != None: # seqres
            #     uniprot = id1
            # else:
            #     region = id1

        elif line.startswith("# 2:"):
            id2 = line.split(":")[1].strip()

            # if regexpUniProt.match(id2) != None: # seqres
            #     uniprot = id2
            # else:
            #     region = id2

        elif line.startswith("# Identity:"):
            identity = line.split(":")[1].strip()
            identity = identity.split('(')[0].strip()
        elif line.startswith("# Similarity:"):
            similarity = line.split(":")[1].strip()
            similarity = similarity.split('(')[0].strip()
        elif line.startswith("# Gaps:"):
            gaps = line.split(":")[1].strip()
            gaps = gaps.split('(')[0].strip()
        elif line.startswith("# Score:"):
            score = line.split(":")[1].strip()


        elif id1 != '' and line.startswith(id1 + ' '):
            sequence = line.split(id1 + ' ')[1]
            seq1 += ''.join([i for i in sequence if not i.isdigit()]).strip()
        elif id2 != '' and line.startswith(id2  + ' '):
            sequence = line.split(id2  + ' ')[1]
            seq2 += ''.join([i for i in sequence if not i.isdigit()]).strip()


        elif region != '' and line.startswith(region.split('_')[0]):

            sequence = line.split('\t')
            if len(sequence) < 2:
                sequence = line.split('    ')[1]
            seq2 += ''.join([i for i in sequence if not i.isdigit()]).strip()

        elif line.startswith("#="):
            count += 1
            if count == 2:
                flag = True
        if flag:
            if line.startswith("#-"):
                flag = False
            alignments += line.replace('\n', '/n')

    if inf is not sys.stdin:
        inf.close()

    len1 = str(len(seq1.replace("-", "")))
    len2 = str(len(seq2.replace("-", "")))


    with open(cfg.data['sequences_regions'] + '/search_in_disprot.tsv', 'r') as input:

        for row in input:

            if id1 == row.split('\t')[0]:

            # if uniprot == row.split('\t')[0]: # seqres
                start = row.split('\t')[6]
                end = row.split('\t')[7]
                disprotid1 = row.split('\t')[4]
                reg1.append(row.split('\t')[5] + '_' + start + '-' + end)

            elif id2 == row.split('\t')[0]:
                start = row.split('\t')[6]
                end = row.split('\t')[7]
                disprotid2 = row.split('\t')[4]
                reg2.append(row.split('\t')[5] + '_' + start + '-' + end)
    reg1 = ','.join(reg1)
    reg2 = ','.join(reg2)


    # out = "\t".join([disprotid1, uniprot, region.split('_')[0], len1, len2, identity, similarity, gaps, score, seq1, seq2, reg1, region]) # seqres
    out = "\t".join([disprotid1, disprotid2, id1, id2, len1, len2, identity, similarity, gaps, score, alignments, reg1, reg2])


    if len(out) > 0:
        print(out)  # script output

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
