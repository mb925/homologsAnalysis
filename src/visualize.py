#!/usr/bin/env python3
import math
import json
import seaborn as sns
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
import config as cfg
from collections import Counter
import statistics


def main():
    # visualize_overlap_identity()
    # visualize_overlap_identity_hist()
    # visualize_alignments()
    # generate_sqv_input('inconsistent_regions')
    # generate_sqv_input('consistent_regions')
    # generate_sqv_input('question_marks')
    # generate_sqv_input('transferable_regions')
    # visualize_pairs(['inconsistent_regions', 'transferable_regions', 'question_marks', 'consistent_regions'])
    visualize_pairs_regions(['inconsistent_regions', 'transferable_regions', 'question_marks', 'consistent_regions'])
    # compare_methods(['transferable_regions', 'inconsistent_regions', 'question_marks', 'consistent_regions'], 'term', 'Disorder function')
    # compare_methods(['transferable_regions', 'inconsistent_regions', 'question_marks', 'consistent_regions'], 'ec', 'Structural state')
    # compare_methods(['transferable_regions', 'inconsistent_regions', 'question_marks', 'consistent_regions'], 'term',
    #                 'Interaction partner')


def compare_methods(dataset_list, evaluated_column, term_namespace):
    df2 = pd.read_csv(cfg.data['sequences_regions'] + '/search_in_disprot.tsv', sep='\t')
    df2 = df2.loc[(df2.term_namespace == term_namespace)]
    if evaluated_column != 'ec':
        df2 = df2.drop_duplicates(subset=['acc', 'start', 'end'], keep=False)
    methods_match_df = pd.DataFrame()
    class_pairs = []
    methods_match = []
    techniques_df = pd.DataFrame(columns=['id1', 'id2', evaluated_column, 'dataset', 'methods_match'])
    for dataset in dataset_list:
        df1 = pd.read_csv(cfg.data['visualizing'] + '/alignment-files/' + dataset + '.csv', sep='\t').groupby(
            ['id1', 'id2']).size()
        reset_df = df1.reset_index()

        for index, row in reset_df.iterrows():
            m1 = df2.loc[(df2.acc == row['id1'])][evaluated_column]
            m2 = df2.loc[(df2.acc == row['id2'])][evaluated_column]
            m1 = m1.to_numpy()
            m2 = m2.to_numpy()

            methods_match_value = len(set(m1).intersection(set(m2)))
            if methods_match_value != 0:
                total = len(set(m1)) + len(set(m2))
                methods_match_value = methods_match_value / total
            methods_match_value = float('%.5f' % methods_match_value)
            methods_match.append(methods_match_value)

            for el1 in m1:
                techniques_df = techniques_df.append(
                    {'id1': row['id1'], 'id2': row['id2'], evaluated_column: el1, 'dataset': dataset.split('_')[0],
                     'methods_match': methods_match_value},
                    ignore_index=True)
            for el2 in m2:
                techniques_df = techniques_df.append(
                    {'id1': row['id1'], 'id2': row['id2'], evaluated_column: el2, 'dataset': dataset.split('_')[0],
                     'methods_match': methods_match_value},
                    ignore_index=True)

            class_pairs.append(dataset.split('_')[0])

        methods_list = techniques_df.loc[(techniques_df.dataset == dataset.split('_')[0])][evaluated_column]
        letter_counts = Counter(methods_list)
        methods_count = pd.DataFrame.from_dict(letter_counts, orient='index')
        methods_count.plot(kind='bar', title=dataset + ' ' + term_namespace)
        plt.tight_layout()
        # plt.savefig(cfg.data['visualizing'] + '/statistics/count/' + evaluated_column + '_' + term_namespace + '_' +
        #             dataset.split('_')[0] + '.png')

    # methods_match_df['matching_' + evaluated_column] = methods_match
    # methods_match_df['dataset'] = class_pairs
    #
    # sns.boxplot(x=methods_match_df["dataset"], y=methods_match_df["matching_" + evaluated_column])
    # plt.tight_layout()
    # plt.savefig(cfg.data['visualizing'] + '/statistics/' + evaluated_column + '_' + term_namespace + '_' + 'matching_methods_class.png')

    # df_heatmap = techniques_df.pivot_table(values='methods_match', index=evaluated_column, columns='dataset', aggfunc=np.mean)
    # sns.heatmap(df_heatmap, annot=True, cmap="BuPu", xticklabels=True, yticklabels=True)
    # plt.tight_layout()
    # plt.savefig(cfg.data['visualizing'] + '/statistics/' + evaluated_column + '_' + term_namespace + '_' + 'methods_heatmap.png')

    print(methods_match_df)

def visualize_pairs_regions(dataset_list):
    length, id, values, cl = [], [], [], []
    data = pd.DataFrame({'length': length, 'id': id, 'value': values, 'class': cl})
    for dataset in dataset_list:
        with open(cfg.data['visualizing'] + '/alignment-files/' + dataset + '.csv', 'r') as file:
            counts1 = {}
            counts2 = {}
            counts_id1 = []
            counts_id2 = []
            pair = ''
            count1 = 0
            count2 = 0
            next(file)
            flag1 = False
            flag2 = False
            for line in file:
                if line.split('\t')[0] + '_' + line.split('\t')[1] != pair:

                    #id1
                    if len(counts_id1) > 0:
                        counts1[line.split('\t')[0]] = statistics.mean(counts_id1)
                    else:
                        counts1[line.split('\t')[0]] = 0
                    count1 = 0
                    flag1 = False
                    counts_id1 = []

                    #id2
                    if len(counts_id2) > 0:
                        counts2[line.split('\t')[1]] = statistics.mean(counts_id2)
                    else:
                        counts2[line.split('\t')[1]] = 0
                    count2 = 0
                    flag2 = False
                    counts_id2 = []

                    pair = line.split('\t')[0] + '_' + line.split('\t')[1]

                else:
                    # id1
                    if line.split('\t')[6] == '1.0':
                        count1 += 1
                        flag1 = True
                    else:
                        if flag1:
                            counts_id1.append(count1)
                            count1 = 0
                            flag1 = False

                    # id2
                    if line.split('\t')[7] == '1.0':
                        count2 += 1
                        flag2 = True
                    else:
                        if flag2:
                            counts_id2.append(count2)
                            count2 = 0
                            flag2 = False

            for el in counts1:
                print(el)
                id.append(el)
                length.append('length_1')
                values.append(counts1[el])
                cl.append(dataset.split('_')[0])

            for el in counts2:
                id.append(el)
                length.append('length_2')
                values.append(counts2[el])
                cl.append(dataset.split('_')[0])

    data['length'] = length
    data['id'] = id
    data['value'] = values
    data['class'] = cl
    sns.violinplot(data=data, y='value', split=True, hue='length', x='class', cut=0)
    plt.show()
    # plt.savefig(cfg.data['visualizing'] + '/statistics/regions_length_class.png')


def visualize_pairs(dataset_list):
    length, id, values, cl = [], [], [], []
    data = pd.DataFrame({'length': length, 'id': id, 'value': values, 'class': cl})
    for dataset in dataset_list:
        df = pd.read_csv(cfg.data['visualizing'] + '/alignment-files/' + dataset + '.csv', sep='\t')
        df_id1 = df.loc[df.s1 != '-']
        df_id2 = df.loc[df.s2 != '-']
        size_df_id1 = df_id1.groupby(['id1', 'id2']).size()
        size_df_id2 = df_id2.groupby(['id1', 'id2']).size()
        id1, id2 = [], []
        for index, value in size_df_id1.items():
           if index[0] not in id1:
               id.append(index[0])
               length.append('length_1')
               id1.append(index[0])
               values.append(value)
               cl.append(dataset.split('_')[0])
        for index, value in size_df_id2.items():
            if index[0] not in id2:
                id.append(index[1])
                length.append('length_2')
                id2.append(index[1])
                values.append(value)
                cl.append(dataset.split('_')[0])

    data['length'] = length
    data['id'] = id
    data['value'] = values
    data['class'] = cl

    sns.violinplot(data=data, y='value', split=True, hue='length', x='class')

    plt.savefig(cfg.data['visualizing'] + '/statistics/protein_length_class.png')


def generate_sqv_input(dataset):
    df = pd.read_csv(cfg.data['visualizing'] + '/alignment-files/' + dataset + '.csv', sep='\t')

    ov_sh = pd.read_csv(cfg.data['visualizing'] + '/overlap-shortest.csv', sep='\t')

    dis = pd.read_csv(cfg.data['sequences_regions'] + '/search-disprot-filtered.tsv', sep='\t')

    ids = df[['id1', 'id2']].drop_duplicates(subset=['id1', 'id2'])

    data = {}
    for el in ids.values:
        example = df.loc[(df.id1 == el[0]) & (df.id2 == el[1])]

        print(example)
        i = 0
        region1 = ''
        sequence1 = ''
        sequence2 = ''
        sequence_middle = ''
        region2 = ''

        while i < len(example['s1'].values):
            if math.isnan(example['d1'].values[i]):
                region1 += str('-')
            else:
                region1 += str(int(example['d1'].values[i]))
            sequence1 += example['s1'].values[i]
            sequence2 += example['s2'].values[i]
            sequence_middle += example['match'].values[i]
            if math.isnan(example['d2'].values[i]):
                region2 += str('-')
            else:
                region2 += str(int(example['d2'].values[i]))
            i += 1

        seqs = []
        row = ov_sh.loc[(ov_sh.id1 == el[0]) & (ov_sh.id2 == el[1])]
        print(el[0])
        # print(el[1])

        organism1 = (dis.loc[dis['acc'] == el[0]])['organism'].iloc[0]
        organism2 = (dis.loc[dis['acc'] == el[1]])['organism'].iloc[0]

        print(organism1)
        print(organism2)

        values = {'overlap': row['overlap-similar'].values[0], 'shortest_region': row['shortest_region'].values[0],
                  'identity': row['identity'].values[0]}
        seqs.append({'sequence': region1, 'id': 1})
        seqs.append({'sequence': sequence1, 'id': 2, 'label': organism1})
        seqs.append({'sequence': sequence_middle, 'id': 3})
        seqs.append({'sequence': sequence2, 'id': 4, 'label': organism2})
        seqs.append({'sequence': region2, 'id': 5})

        data[el[0] + '-' + el[1]] = {'data': seqs, 'values': values}

    with open(cfg.data['visualizing'] + '/alignment-files/sqv/' + dataset + '.json', 'w') as json_file:
        json.dump(data, json_file)


def visualize_overlap_identity():
    ident, region_identity, region_overlap, union_df, local_region_overlap, local_region_identity = [], [], [], [], [], []


    df = pd.read_csv(cfg.data['rearrange'] + '/all-dataframe.tsv', sep='\t')
    identity_df = pd.read_csv(cfg.data['clustering'] + '/components-filtered.tsv', sep='\t')

    ids_identity = identity_df[['id1', 'id2', 'identity']].sort_values(by=['id1'])
    union_df = df.loc[(df.d2 == 1) | (df.d1 == 1)].groupby(['id1', 'id2']).size().to_frame('union')
    union_region1 = df.loc[(df.d1 == 1)].groupby(['id1', 'id2']).size().to_frame('union-region1').sort_values(
        by=['id1', 'id2'])

    union_region2 = df.loc[(df.d2 == 1)].groupby(['id1', 'id2']).size().to_frame('union-region2').sort_values(
        by=['id1', 'id2'])
    min_lens = df.groupby(['id1', 'id2']).size()

    join = union_region2.merge(union_region1, on=['id1', 'id2'], how='right')
    # print(union_region1)
    # print(union_region2)

    m1 = pd.merge(ids_identity, union_df, on=['id1', 'id2'])
    # print(m1)
    overlap_similar = df.loc[(df.d2 == 1) & (df.d1 == 1)].groupby(['id1', 'id2']).size().to_frame('overlap-similar')
    m2 = pd.merge(m1, overlap_similar, on=['id1', 'id2'], how='outer').sort_values(by=['id1', 'id2'])
    # print(m2)
    overlap_match = df.loc[(df.d2 == 1) & (df.d1 == 1) & (df.s1 == df.s2)].groupby(['id1', 'id2']).size().to_frame(
        'overlap-match')
    m3 = pd.merge(m2, overlap_match, on=['id1', 'id2'], how='outer')
    m4 = pd.merge(m3, union_region1, on=['id1', 'id2'], how='outer')
    m5 = pd.merge(m4, union_region2, on=['id1', 'id2'], how='outer')

    for value in m5['identity'].values:
        num = int(value.split('/')[0])
        union = int(value.split('/')[1])

        value = 0
        if num != 0:
            value = (num / union) * 100
        ident.append(round(value, 2))

    m5['identity'] = m5['identity'].replace(m5['identity'].values, ident)

    i = 0

    while i < len(m5['overlap-similar'].values):
        ov_sim = 0
        ov_match = 0
        num_sim = 0
        num_match = 0
        union = m5['union'][i]

        if not np.isnan(m5['overlap-similar'][i]):
            num_sim = int(m5['overlap-similar'][i])
        if not np.isnan(m5['overlap-match'][i]):
            num_match = int(m5['overlap-match'][i])

        if num_sim != 0:
            ov_sim = round((num_sim / union), 2)

        if num_match != 0:
            ov_match = round((num_match / union), 2)

        region_overlap.append(ov_sim)
        region_identity.append(ov_match)

        if num_sim != 0:
            # local region overlap: matches and mismatches over the length of the shortest region

            minimum = min(m5['union-region1'][i], m5['union-region2'][i])
            print(minimum, ' ', union)
            if minimum == union:
                print('hey!', ' ', minimum, ' ', union, ' ', num_sim)
                print(m5.iloc[i])
            local_region_overlap.append(round(num_sim / minimum, 2))
        else:
            local_region_overlap.append(0)

        if num_match != 0:
            # region overlap / shortest region
            minimum = min(m5['union-region1'][i], m5['union-region2'][i])
            local_region_identity.append(round(num_match / minimum, 2))
        else:
            local_region_identity.append(0)

        # print(round(num_sim / min(m5['union-region1'][i], m5['union-region1'][i]), 2))
        i += 1

    m5['overlap-similar/union'] = region_overlap
    # print(m5['union-region2'].values)
    # print(m5['union-region1'].values)
    # print(m5['overlap-similar'])
    # print(shortest_region)
    m5.to_csv(cfg.data['visualizing'] + '/overlap-union-identity.csv', index=None, sep='\t')
    ov_sh = pd.DataFrame({'id1': m5['id1'], 'id2': m5['id2'], 'identity': m5['identity'],
                          'local_region_overlap': local_region_overlap,
                          'region_overlap': region_overlap})

    # print(m5)

    ov_sh.to_csv(cfg.data['visualizing'] + '/overlap-shortest.csv', index=None, sep='\t')

    fig = plt.figure(figsize=(35, 10), dpi=50)

    fig.suptitle('overlap/union', fontsize=24)

    ax1 = fig.add_subplot(1, 3, 1)
    ax1.set_title('colored by region overlap')
    ax1.set_ylabel('region overlap', fontsize=20)
    ax1.set_xlabel('global sequence identity', fontsize=20)
    cm = plt.cm.get_cmap('seismic')
    sc = ax1.scatter(ident, region_overlap, c=m5['overlap-similar'], cmap=cm)
    plt.colorbar(sc)

    # ax2 = fig.add_subplot(1, 3, 2)
    # sc2 = ax2.scatter(ident, region_overlap, c=local_region_overlap, cmap=cm)
    # plt.colorbar(sc2)
    #
    # ax2.set_xlabel('global sequence identity', fontsize=20)
    # ax2.set_ylabel('region overlap', fontsize=20)
    # ax2.set_title('colored by local region overlap')
    #
    # ax3 = fig.add_subplot(1, 3, 3)
    # sc3 = ax3.scatter(ident, region_overlap, c=local_region_identity, cmap=cm)
    # plt.colorbar(sc3)
    #
    # ax3.set_xlabel('global sequence identity', fontsize=20)
    # ax3.set_ylabel('region overlap', fontsize=20)
    # ax3.set_title('colored by local region identity')
    plt.show()

    # fig.savefig(cfg.data['visualizing'] + '/global_identity-region_overlap.png')
    #
    # cm = plt.cm.get_cmap('seismic')
    # fig = plt.figure(figsize=(20, 10), dpi=50)
    #
    # # fig.suptitle('region identity and region overlap shortest', fontsize=24)
    # ax1 = fig.add_subplot(1, 2, 1)
    # ax1.set_xlabel('region overlap', fontsize=20)
    # ax1.set_ylabel('region identity', fontsize=20)
    # ax1.scatter(over_mis, over_mat)
    #
    # ax2 = fig.add_subplot(1, 2, 2)
    # ax2.set_xlabel('region overlap', fontsize=20)
    # ax2.set_ylabel('local region identity', fontsize=20)
    # ax2.scatter(over_mis, shortest_region)
    #
    # plt.show()
    # fig.savefig(cfg.data['visualizing'] + '/region-identity_local-region-identity.png')


def normalize_data(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))


def visualize_overlap_identity_hist():
    data = []
    dict = {}
    ids = []
    data_overlap = []

    df = pd.read_csv(cfg.data['visualizing'] + '/overlap-union-identity.csv', sep='\t')
    for i, id1 in enumerate(df['id1']):
        data_overlap.append(df['overlap-similar/union'][i])

        if id1 in dict:
            dict[id1].append(df['identity'][i])

        else:
            dict[id1.split('_')[0]] = []
            dict[id1.split('_')[0]].append(df['identity'][i])
    for i, id2 in enumerate(df['id2']):
        if id2 in dict:
            dict[id2].append(df['identity'][i])

        else:
            dict[id2.split('_')[0]] = []
            dict[id2.split('_')[0]].append(df['identity'][i])
    print(dict)
    for key in dict:
        max_val = max(dict[key])
        ids.append(key)
        data.append(round(max_val, 2))

    print(sorted(data_overlap))

    # plt.hist(data, bins=20) # !!!do not run both plots together or the results will overlap
    # plt.xlabel('max global sequence identity')
    # plt.ylabel('# proteins')
    # # plt.show()
    # plt.gcf().savefig(cfg.data['visualizing'] + '/max_global_sequence_identity.png')

    plt.xlabel('region overlap')
    plt.ylabel('# pairs')

    plt.hist(data_overlap, bins=20)
    plt.gcf().savefig(cfg.data['visualizing'] + '/region_overlap_distribution.png')


def visualize_alignments():
    df = pd.read_csv(cfg.data['rearrange'] + '/all-dataframe.tsv', sep='\t')
    m5 = pd.read_csv(cfg.data['visualizing'] + '/overlap-union-identity.csv', sep='\t')
    inconsistent_regions = \
    m5.loc[(m5['overlap-similar/union'] <= 0.5) & (m5['overlap-similar/union'] >= 0.1) & (m5['identity'] >= 50)][
        ['id1', 'id2']]
    question_marks = m5.loc[(m5['overlap-similar/union'] <= 0.2) & (m5['identity'] >= 50)][['id1', 'id2']]
    transferable_regions = m5.loc[(m5['overlap-similar/union'] >= 0.5) & (m5['identity'] <= 50)][['id1', 'id2']]
    consistent_regions = m5.loc[(m5['overlap-similar/union'] >= 0.5) & (m5['identity'] >= 50)][['id1', 'id2']]
    print(m5)
    inconsistent_dataset = df.merge(inconsistent_regions, on=['id1', 'id2'], how='right')
    inconsistent_dataset.to_csv(cfg.data['visualizing'] + '/alignment-files/inconsistent_regions.csv', index=None,
                                sep='\t')
    question_marks_dataset = df.merge(question_marks, on=['id1', 'id2'], how='right')
    question_marks_dataset.to_csv(cfg.data['visualizing'] + '/alignment-files/question_marks.csv', index=None, sep='\t')
    transferable_regions_dataset = df.merge(transferable_regions, on=['id1', 'id2'], how='right')
    transferable_regions_dataset.to_csv(cfg.data['visualizing'] + '/alignment-files/transferable_regions.csv',
                                        index=None, sep='\t')
    consistent_regions_dataset = df.merge(consistent_regions, on=['id1', 'id2'], how='right')
    consistent_regions_dataset.to_csv(cfg.data['visualizing'] + '/alignment-files/consistent_regions.csv', index=None,
                                      sep='\t')


if __name__ == '__main__':
    sys.exit(main())
