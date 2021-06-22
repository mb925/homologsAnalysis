import os
import pathlib
import config as cfg

curpath = cfg.fasta['sequences']
curpathreg = cfg.fasta['regions']
outfile = cfg.alignments['global'] + '/global_triangular_matrix.csv'
# outfilereg =  cfg.alignments['global'] + '/regseq_triangular_matrix.csv'


# outfile = 'alignments_global_needle/regseq_triangular_matrix.csv'
# create the couples
def create_pairwise_regseq(dir1, dir2, outfile):

    with open(outfile, 'w+') as regseq_triangular_matrix:
        for uniprot in os.listdir(dir1):
            for region in os.listdir(dir2):
                regseq_triangular_matrix.write(uniprot.split('.')[0] + '\t' + region.split('.')[0] + '\n')



# outfile = 'alignments_global_needle/global_triangular_matrix.csv'
# create the couples
# creating a triangular matrix
def create_pairwise_sequences(dir1, outfile): # TODO: remove diagonal from matrix!
    pairwise_sequences = {}

    with open(outfile, 'w') as global_triangular_matrix:

        for filename1 in os.listdir(dir1):

            i = os.listdir(curpath).index(filename1)
            # print(filename1)
            # for filename2 in os.listdir(curpath)[i:-1]:
                # print(filename2)

            for filename2 in os.listdir(dir1): # for squared matrix
                # if pairwise_sequences.get(filename1.split(".")[0], -1) == -1:
                if filename1.split(".")[0] != filename2.split(".")[0]:
                    if pairwise_sequences.get(filename1.split(".")[0] + '_' + filename2.split(".")[0], -1) == -1:
                        if pairwise_sequences.get(filename2.split(".")[0] + '_' + filename1.split(".")[0], -1) == -1:
                            pairwise_sequences[filename1.split(".")[0] + '_' + filename2.split(".")[0]] = []
                            global_triangular_matrix.write(filename1.split(".")[0] + '\t' + filename2.split(".")[0] + '\n')

                # else:
                #     if filename1.split(".")[0] != filename2.split(".")[0]:
                #         pairwise_sequences[filename1.split(".")[0]].append(filename2.split(".")[0])
                        # global_triangular_matrix.write(filename1.split(".")[0] + '\t' + filename2.split(".")[0] + '\n')


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # create_pairwise_sequences(curpath, outfile)
    # create_pairwise_regseq(curpath, curpathreg, outfilereg)
    create_pairwise_sequences(curpath, outfile)
