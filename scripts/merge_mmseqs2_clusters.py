from Bio import SeqIO
import argparse
from pathlib import Path
from natsort import natsorted
import pickle
import os

def merge_mmseqs2(indir, clusters):
    rep_files = []
    for path in Path(indir).glob("*_cluster.tsv"):
        # Print the path (file or directory) to the console
        rep_files.append(path)

    # sort based on batch number
    rep_files = natsorted(rep_files)

    # fastas are hierarchically clustered, so need to extract the representatives in each file
    # iteratively and then determine which cluster they form in the next file
    rep_to_cluster = {}
    cluster_to_rep = {}
    for rep_file in rep_files:
        with open(rep_file, "r") as f:
            while True:
                line = f.readline()
                if not line:
                    break

                split_line = line.rstrip().split("\t")
                rep = split_line[0]
                seq = split_line[1]

                # if rep has not been seen before, add to new cluster
                if rep not in rep_to_cluster:
                    rep_to_cluster[rep] = rep
                    cluster_to_rep[rep] = set()
                    cluster_to_rep[rep].add(rep)

                # if sequence is in rep_to_cluster, means it has been 
                # clustered with a new represenative
                if seq in rep_to_cluster and seq != rep:
                    # update old reps
                    reps_set = cluster_to_rep[seq]

                    for prev_rep in reps_set:
                        rep_to_cluster[prev_rep] = rep

                    # account for duplicated IDs
                    try:
                        cluster_to_rep[rep].update(reps_set)
                        del cluster_to_rep[seq]
                    except KeyError:
                        pass
                # elif seq != rep:
                #     # add sequence to cluster
                #     cluster_to_rep[rep].add(seq)
        print("Finished: {}".format(rep_file))

    print("No. clusters: {}".format(str(len(cluster_to_rep))))
    output_dicts = (rep_to_cluster, cluster_to_rep)

    with open(clusters, 'wb') as handle:
        pickle.dump(output_dicts, handle, protocol=pickle.HIGHEST_PROTOCOL)

merge_mmseqs2_clusters(snakemake.input.indir, snakemake.output.clusters)