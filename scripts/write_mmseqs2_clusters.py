from Bio import SeqIO
import argparse
from pathlib import Path
from natsort import natsorted
import pickle
import os

def write_mmseqs2_clusters(rep_files, outfile, clusters):
    # sort based on batch number
    rep_files = natsorted(rep_files)

    # read in cluster file
    with open(clusters, 'rb') as handle:
        output_dicts = pickle.load(handle)
        rep_to_cluster, cluster_to_rep = output_dicts
    
    # write output
    with open(outfile, "w") as o:
        for rep_file in rep_files:
            current_rep = None
            with open(rep_file, "r") as f:
                while True:
                    line = f.readline()
                    if not line:
                        break
                    split_line = line.rstrip().split("\t")

                    rep = split_line[0]
                    seq = split_line[1]

                    new_rep = rep_to_cluster[rep]

                    o.write(new_rep + "\t" + seq + "\n")

            print("Finished: {}".format(rep_file))

write_mmseqs2_clusters(snakemake.input.infiles, snakemake.output.outfile, snakemake.input.clusters)