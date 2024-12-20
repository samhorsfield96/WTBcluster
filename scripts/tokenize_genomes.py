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

def get_gff(directory):
    files = (str(p.resolve()) for p in Path(directory).rglob("*.gff*"))
    return list(files)

def chunks(l, n):
    """Yield n number of striped chunks from l."""
    for i in range(0, n):
        yield l[i::n]

def tokenise_gff(files_list, outfile, gene_tokens):
    with open(outfile, "w") as o:
        for gff in files_list:
            basename = os.path.basename(gff)
            with open(gff, "r") as f:
                tokenised_genome = []
                current_contig = None
                while True:
                    line = f.readline()
                    if not line or line == "##FASTA\n":
                        break
                    
                    # skip commented lines
                    if line[0] == "#":
                        continue

                    split_line = line.rstrip().split("\t")
                    type = split_line[2]
                    # if type == "region":
                    #     # add space between contigs as synteny is unknown
                    #     if len(tokenised_genome) > 0:
                    #         tokenised_genome.append("_")
                    # else:
                    
                    gene_strand = True if split_line[6] == "+" else False
                    split_gene_id = split_line[-1].split(";")[0].replace("ID=", "")
                    
                    contig_ID = split_gene_id[0].zfill(5)

                    # add contig end
                    if contig_ID != current_contig:
                        if len(tokenised_genome) > 0:
                            tokenised_genome.append("_")
                        current_contig = contig_ID

                    # build gene id to search in dictionary
                    gene_name = split_gene_id

                    gene_token = gene_tokens.get(gene_name.encode())
                    if gene_token is not None:
                        gene_token = gene_token.decode()
                        # multiply by minus 1 for negative strand
                        if not gene_strand:
                            gene_token = "-" + gene_token
                        
                        tokenised_genome.append(str(gene_token))
            
            tokenised_genome_str = " ".join(tokenised_genome)
            o.write(basename + "\t" + tokenised_genome_str + "\n")

def tokenize_genomes(batch_file, outfile, gene_tokens_db):

    opts = rocksdb.Options()
    opts.max_open_files = 300000000
    opts.max_bytes_for_level_base = 209715200 #(default = 10485760)
    opts.target_file_size_base = math.ceil(opts.max_bytes_for_level_base / 10) #(default = 2097152)
    opts.target_file_size_multiplier = 2 #(default = 1)

    gene_tokens = rocksdb.DB(gene_tokens_db, opts, read_only=True)

    # generate list of integers to represent genome
    files_list = []

    with open(batch_file, "r") as o:
        while True:
            line = o.readline()
            if not line:
                break
            files_list.append(line.rstrip())

    file_index = batch_file.split("mmseqs2_batch_N")[-1].replace("_gff.txt", "")

    tokenise_gff(files_list, outfile, gene_tokens)

        
tokenize_genomes(snakemake.input.batch_file, snakemake.output.outfile, snakemake.input.out_db)