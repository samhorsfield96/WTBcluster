from Bio import SeqIO
import os
import subprocess
import argparse
import sys
from collections import defaultdict
import pathlib

def get_options():
    description = "Query sequence in WTBclusters MMseqs2 output"
    parser = argparse.ArgumentParser(description=description,
                                        prog='python query_mmseqs2.py')
    IO = parser.add_argument_group('Input/options.out')
    IO.add_argument('--query',
                    required=True,
                    help='FASTA File containing sequences search for.')
    IO.add_argument('--reps',
                    required=True,
                    help='FASTA File merged representatives from MMseqs2.')
    IO.add_argument('--mmseqs2-params',
                    default="",
                    help='String of MMseqs2 parameters for mmseqs easy-search in quotation ("") marks. Default runs default params.'
                        '\nExample: "--min-seq-id 0.5 -c 0.5 --seq-id-mode 0 --threads 4 --cov-mode 0 --alignment-mode 3"')
    IO.add_argument('--mmseqs2-merged',
                    required=True,
                    help='Merged MMseqs2 all_clusters.tsv.')
    IO.add_argument('--mmseqs2-all-fasta',
                    required=True,
                    help='Directory containing MMseqs all_seqs.fasta files.')
    IO.add_argument('--tmp',
                    required=True,
                    help='Temporary directory.')
    IO.add_argument('--outpref',
                    required=True,
                    help='Output preferences')

    return parser.parse_args()

def main():
    options = get_options()
    query = options.query
    reps = options.reps
    mmseqs2_params = options.mmseqs2_params.replace('"', '')
    #print(mmseqs2_params)
    mmseqs2_merged = options.mmseqs2_merged
    mmseqs2_all_fasta = options.mmseqs2_all_fasta
    outpref = options.outpref
    tmp = options.tmp

    # run initial easy-search
    search_output = outpref + ".m8"
    try:
        subprocess.run([
            "mmseqs", "easy-search", query, reps, search_output, tmp, *mmseqs2_params.split(" ")
        ], check=True)

        print(f"Clustering results saved to {search_output}")

    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}")
        sys.exit(1)
    
    # read in output file
    print(f"Reading alignments {search_output}...")
    centroid_dict = defaultdict(set)
    with open(search_output, "r") as f:
        while True:
            line = f.readline()
            if not line:
                break

            split_line = line.rstrip().split("\t")

            query_id = split_line[0]
            centroid_id = split_line[1]
            centroid_dict[centroid_id].add(query_id)
    #print("centroid_dict")
    #print(centroid_dict)

    # search for hit sequences in merged mmseqs files
    cluster_dict = defaultdict(set)
    seq_to_centroid_dict = {}
    print(f"Reading all clusters {mmseqs2_merged}...")
    with open(mmseqs2_merged, "r") as f:
        while True:
            line = f.readline()
            if not line:
                break

            split_line = line.rstrip().split("\t")
            centroid_id = split_line[0]
            seq_id = split_line[1]

            # add hit sequence to final output dict
            if centroid_id in centroid_dict:
                seq_to_centroid_dict[seq_id] = centroid_id
                for query_id in centroid_dict[centroid_id]:
                    cluster_dict[seq_id].add(query_id)
    #print("cluster_dict")
    #print(cluster_dict)
    #print("seq_to_centroid_dict")
    #print(seq_to_centroid_dict)

    # go through all files in mmseqs2_all_fasta directory and find matches
    mmseqs2_all_fasta_dir = pathlib.Path(mmseqs2_all_fasta)
    mmseqs2_all_fasta_list = list(mmseqs2_all_fasta_dir.glob("*_all_seqs.fasta"))

    # search for hits in all fasta files
    final_dict = defaultdict(set)
    for input_file in mmseqs2_all_fasta_list:
        print(f"Reading FASTA {input_file}...")
        fasta_sequences = SeqIO.parse(open(input_file),'fasta')
        for fasta in fasta_sequences:
            name, description, sequence = fasta.id, fasta.description, str(fasta.seq)

            # skip empty fields
            if len(sequence) == 0:
                continue

            if name in cluster_dict:
                for query_id in cluster_dict[name]:
                    final_dict[query_id].add((name, description, sequence))
    
    #print("final_dict")
    #print(final_dict)

    # write the fasta files
    print(f"Printing output files...")
    for query_id, hit_list in final_dict.items():
        with open(outpref + "_" + query_id + ".fasta", "w") as o:
            for name, description, sequence in hit_list:
                centroid_id = seq_to_centroid_dict[name]
                o.write(">" + description + " cluster: " + centroid_id + "\n" + sequence + "\n")
    
    sys.exit(0)

                
















if __name__ == "__main__":
    main()