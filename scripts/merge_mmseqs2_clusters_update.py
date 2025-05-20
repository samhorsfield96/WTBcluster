from Bio import SeqIO
import argparse
from pathlib import Path
from natsort import natsorted
import pickle
import os

def merge_mmseqs2(rep_files, prev_clusters, clusters):
    # read in cluster file
    with open(prev_clusters, 'rb') as handle:
        output_dicts = pickle.load(handle)
        prev_rep_to_cluster, prev_cluster_to_rep = output_dicts

    # sort based on batch number
    rep_files = natsorted(rep_files)


    rep_to_cluster = {}
    cluster_to_rep = {}

    merged_old_reps = set()
    
    # keep track of any representatives that have swapped
    rep_swap = {}

    # iterate over each rep_file, creating tmp file for each and overwritting
    # to ensure that reps in original clustering as kept as reps this time
    for rep_file in rep_files:
        with open(rep_file, "r") as f:
            while True:
                line = f.readline()
                if not line:
                    break

                split_line = line.rstrip().split("\t")
                rep = split_line[0]
                seq = split_line[1]

                key = None
                entry = None

                # fine, old rep matches itself
                if rep in prev_rep_to_cluster and seq in prev_rep_to_cluster and rep == seq:
                    key, entry = (rep, seq)
                
                # fine, new sequence added to old rep
                elif rep in prev_rep_to_cluster and seq not in prev_rep_to_cluster:
                    key, entry = (rep, seq)
                
                # fine, completely new cluster
                elif rep not in prev_rep_to_cluster and seq not in prev_rep_to_cluster:
                    key, entry = (rep, seq)

                # problem, merging of two old reps
                elif rep in prev_rep_to_cluster and seq in prev_rep_to_cluster and rep != seq:
                    merged_old_reps.add((rep, seq))
                    continue

                # problem, new seq added as new rep, potential merge of multiple old reps
                elif rep not in prev_rep_to_cluster and seq in prev_rep_to_cluster:
                    # if already present, indicates multiple old reps have been merged
                    if seq not in rep_swap:
                        rep_swap[rep] = seq
                    else:
                        # add "lowest" to be the representative
                        if rep_swap[rep] < seq:
                            rep_swap[rep] = seq
                        
                        print(f"MERGE new_rep: {rep} old_rep {seq}")
                    
                    key, entry = (rep_swap[rep], rep)
                    
                # if rep has not been seen before, add to new cluster
                if key not in rep_to_cluster:
                    rep_to_cluster[key] = key
                    cluster_to_rep[key] = set()
                    cluster_to_rep[key].add(key)

                # if sequence is in cluster_to_rep, means it has been 
                # clustered with a new represenative
                if entry in cluster_to_rep and key != entry:
                    # update old reps
                    reps_set = cluster_to_rep[entry]

                    for prev_rep in reps_set:
                        rep_to_cluster[prev_rep] = key

                    # account for duplicated IDs
                    try:
                        cluster_to_rep[key].update(reps_set)
                        del cluster_to_rep[entry]
                    except KeyError:
                        pass
        print("Finished: {}".format(rep_file))

    print("No. clusters: {}".format(str(len(cluster_to_rep))))
    output_dicts = (rep_to_cluster, cluster_to_rep)

    with open(clusters, 'wb') as handle:
        pickle.dump(output_dicts, handle, protocol=pickle.HIGHEST_PROTOCOL)

merge_mmseqs2(snakemake.input.infiles, snakemake.input.prev_clusters, snakemake.output.clusters)