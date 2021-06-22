import config as cfg
# Press the green button in the gutter to run the script.
import filter
import cluster
import rearrange
import visualize

if __name__ == '__main__':
    # filter >= 30% identity
    filter.filter_identity(30)
    # generate clusters
    cluster.cluster_file()  # useful only when the filtering is large, otherwise data are all connected in a big cluster
    # keep uniprots with a structural state region associated
    filter.filter_structural()
    # filter uniprots from clusters
    cluster.filter_cluster()
    # rearrange data into dataframes
    rearrange.create_tables()
    # calculate union overlap and visualize a graph
    visualize.visualize_overlap_identity()
    # visualize a subset of alignments
    ## study of cases with 80 - 100 % identity
    visualize.visualize_alignments()
