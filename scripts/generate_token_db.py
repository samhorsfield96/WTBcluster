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

def generate_gene_id(unsplit_id):
    name = unsplit_id.split("SAM")[1].split("_")[0].split(".")[0]
    split_gene = unsplit_id.split("_")
    gene_id = split_gene[-1]
    contig_id = split_gene[-2][-5:]

    gene_name = name + "_" + contig_id + "_" + gene_id

    return gene_name

def generate_token_db(cluster_file, out_db, out_reps):

    opts = rocksdb.Options()
    opts.max_open_files = 300000000
    opts.max_bytes_for_level_base = 209715200 #(default = 10485760)
    opts.target_file_size_base = math.ceil(opts.max_bytes_for_level_base / 10) #(default = 2097152)
    opts.target_file_size_multiplier = 2 #(default = 1)

    # dictionary of representative sequences and their token
    reps_dict = {}
    rep_to_token = {}

    # dictionary mapping each gene to a given cluster token
    gene_tokens_db = out_db
    opts.create_if_missing = True
    gene_tokens = rocksdb.DB(gene_tokens_db, opts)

    #start at 0 as cannot assign negative 0
    token = 0
    current_rep = None
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

            rep = generate_gene_id(split_rep)
            seq = generate_gene_id(split_seq)

            # allows use of non-sorted list
            if rep not in rep_to_token:
                token += 1
                reps_dict[token] = rep
                rep_to_token[rep] = token
            
            current_token = rep_to_token[rep]
            
            # add sequence to cluster
            gene_tokens.put(seq.encode(), str(current_token).encode())
            counter += 1
            if counter % 10000000 == 0:
                print("At index: {}".format(counter))

    # save data as pickle
    print("Saving token dictionaries...")

    with open(out_reps, "wb") as f:
        pickle.dump((reps_dict, rep_to_token), f)
    
    print("Saved token dictionaries.")

generate_token_db(snakemake.input.clusters, snakemake.output.out_db, snakemake.output.out_reps)