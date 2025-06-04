import Bio.SeqIO
import pyrodigal
import gzip
import os

from pathlib import Path
import argparse
import os
import pickle
#from multiprocessing import Pool, Manager
#from functools import partial
import rocksdb
import math

def update_token_db(cluster_file, out_db, in_reps, out_reps):

    # read in representative sequences and their token
    with open(in_reps, 'rb') as handle:
        reps_dict, rep_to_token = pickle.load(handle)

    # dictionary mapping each gene to a given cluster token, update in place
    opts = rocksdb.Options()
    opts.max_open_files = 300000000
    opts.max_bytes_for_level_base = 209715200 #(default = 10485760)
    opts.target_file_size_base = math.ceil(opts.max_bytes_for_level_base / 10) #(default = 2097152)
    opts.target_file_size_multiplier = 2 #(default = 1)

    gene_tokens = rocksdb.DB(out_db, opts)

    #start at len(rep_to_token), will produce new token when new one discovered
    token = len(rep_to_token)
    print("Generating token dictionaries...")
    counter = 0
    with open(cluster_file, "r") as f:
        while True:
            line = f.readline()
            if not line:
                break

            split_line = line.rstrip().split("\t")
            split_rep = split_line[0]
            split_seq = split_line[1]

            rep = split_rep
            seq = split_seq

            # if previously clustered rep is the seq, ignore as don't want to reassign
            if seq in rep_to_token:
                continue

            # allows use of non-sorted list
            if rep not in rep_to_token:
                token += 1
                rep_to_token[rep] = token
            
            current_token = rep_to_token[rep]
            reps_dict[current_token] = rep
            
            # add sequence to cluster
            gene_tokens.put(seq.encode(), str(current_token).encode())
            counter += 1
            if counter % 1000000 == 0:
                print("At index: {}".format(counter))

    # save data as pickle
    print("Saving token dictionaries...")

    with open(out_reps, "wb") as f:
        pickle.dump((reps_dict, rep_to_token), f)
    
    print("Saved token dictionaries.")

update_token_db(snakemake.input.clusters, snakemake.input.out_db, snakemake.input.in_reps, snakemake.output.out_reps)