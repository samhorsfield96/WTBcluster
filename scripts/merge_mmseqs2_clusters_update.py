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

    #prev_rep_to_cluster = {"old1": "old1", "old2": "old2", "old3": "old3", "old4": "old4", "old5": "old5", "old6": "old6"}

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
                    #print(f"Option1: {(key, entry)}")
                
                # fine, new sequence added to old rep
                elif rep in prev_rep_to_cluster and seq not in prev_rep_to_cluster:
                    key, entry = (rep, seq)
                    #print(f"Option2: {(key, entry)}")
                
                # fine, completely new cluster
                elif rep not in prev_rep_to_cluster and seq not in prev_rep_to_cluster:
                    key, entry = (rep, seq)
                    #print(f"Option3: {(key, entry)}")

                # problem, merging of two old reps
                elif rep in prev_rep_to_cluster and seq in prev_rep_to_cluster and rep != seq:
                    #merged_old_reps.add((rep, seq))
                    #print(f"Option4: {(rep, seq)}")
                    continue

                # problem, new seq added as new rep, potential merge of multiple old reps
                elif rep not in prev_rep_to_cluster and seq in prev_rep_to_cluster:
                    # if already present, indicates multiple old reps have been merged
                    if rep not in rep_swap:
                        rep_swap[rep] = seq
                    else:
                        # add "lowest" to be the representative
                        old_swap = rep_swap[rep]
                        if rep_swap[rep] < seq:
                            rep_swap[rep] = seq
                        
                        print(f"MERGE new_rep: {rep}, old_rep: {seq}, prev_old_rep: {old_swap}")
                    
                    key, entry = (rep_swap[rep], rep)
                    #print(f"Option5: {(key, entry)}")
                    
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
        #print(f"rep_to_cluster: {rep_to_cluster}")
        #print(f"cluster_to_rep: {cluster_to_rep}")
        print("Finished: {}".format(rep_file))

    print("No. clusters: {}".format(str(len(cluster_to_rep))))
    output_dicts = (rep_to_cluster, cluster_to_rep)

    with open(clusters, 'wb') as handle:
        pickle.dump(output_dicts, handle, protocol=pickle.HIGHEST_PROTOCOL)

#merge_mmseqs2(["/hps/software/users/jlees/shorsfield/software/WTBcluster/scripts/test_batch1.tsv", "/hps/software/users/jlees/shorsfield/software/WTBcluster/scripts/test_batch2.tsv"], None, "test_clustering.txt")
#print("Finished forward")
#merge_mmseqs2(["/hps/software/users/jlees/shorsfield/software/WTBcluster/scripts/test_batch2.tsv", "/hps/software/users/jlees/shorsfield/software/WTBcluster/scripts/test_batch1.tsv"], None, "test_clustering.txt")
#print("Finished reverse")
merge_mmseqs2(snakemake.input.infiles, snakemake.input.prev_clusters, snakemake.output.clusters)