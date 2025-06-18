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
    IO.add_argument('--mmseqs2-params',
                    default="",
                    help='String of MMseqs2 parameters for mmseqs easy-search in quotation ("") marks. Default runs default params.'
                        '\nExample: "--min-seq-id 0.5 -c 0.5 --seq-id-mode 0 --threads 4 --cov-mode 0 --alignment-mode 3"')
    IO.add_argument('--fasta-dir',
                    required=True,
                    help='Directory containing all .fasta files.')
    IO.add_argument('--fasta-pref',
                    default="",
                    help='Optional file prefix for .fasta files.')
    IO.add_argument('--m8_file',
                    default=None,
                    help='Previous results from query_mmseqs2.py.')
    IO.add_argument('--length-discrepancy',
                    default="0.0,0.0",
                    help='Proportional length difference between query and returned proteins. Pass as comma separated list e.g. 0.9,1.1 will return proteins that are 10%/ smaller and larger than the query. If unspecifed, returns all.'
                        'To set just one boundary, set the other to 0.0 e.g. 0.9,0.0 checks only the lower length boundary.')
    IO.add_argument('--no-partials',
                    default="False",
                    action="store_true",
                    help="Don't return partial sequences (i.e. without start and stop codon)")
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
    mmseqs2_params = options.mmseqs2_params.replace('"', '')
    fasta_dir = options.fasta_dir
    outpref = options.outpref
    tmp = options.tmp
    length_discrepancy = options.length_discrepancy
    no_partials = options.no_partials
    m8_file = options.m8_file

    length_discrepancy = [float(x) for x in length_discrepancy.split(",")]
    assert(len(length_discrepancy) == 2)

    # get query sequences and store
    query_dict = {}
    fasta_sequences = SeqIO.parse(open(query),'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        query_dict[name] = sequence

      # go through all files in mmseqs2_all_fasta directory and find matches
    mmseqs2_all_fasta_dir = pathlib.Path(fasta_dir)
    mmseqs2_all_fasta_list = list(mmseqs2_all_fasta_dir.glob("*" + options.fasta_pref + ".fasta"))

    # search for hits in all fasta files
    final_dict = defaultdict(set)
    for input_file in mmseqs2_all_fasta_list:
        # run initial easy-search

        if m8_file == None:
            search_output = outpref + ".m8"
            print(f"Searching FASTA {input_file}...")
            try:
                subprocess.run([
                    "mmseqs", "easy-search", query, input_file, search_output, tmp, *mmseqs2_params.split(" ")
                ], check=True)

                print(f"Clustering results saved to {search_output}")

            except subprocess.CalledProcessError as e:
                print(f"An error occurred: {e}")
                sys.exit(1)

        else: 
            search_output = m8_file
            
        # read in output file
        print(f"Reading alignments {search_output}...")
        cluster_dict = defaultdict(set)
        with open(search_output, "r") as f:
            while True:
                line = f.readline()
                if not line:
                    break

                split_line = line.rstrip().split("\t")

                query_id = split_line[0]
                centroid_id = split_line[1]
                cluster_dict[centroid_id].add(query_id)

        print(f"Reading FASTA {input_file}...")
        fasta_sequences = SeqIO.parse(open(input_file),'fasta')
        for fasta in fasta_sequences:
            name, description, sequence = fasta.id, fasta.description, str(fasta.seq)

            # skip empty fields
            if len(sequence) == 0:
                continue

            # Check if partials need to be returned
            if no_partials and ";partial=00;" not in description:
                continue

            if name in cluster_dict:
                for query_id in cluster_dict[name]:
                    
                    len_query = len(query_dict[query_id])

                    # determine if length of sequence is too short
                    if length_discrepancy[0] > 0.0:
                        if len(sequence) < len_query * length_discrepancy[0]:
                            continue
                    
                    # determine if length of sequence is too long
                    if length_discrepancy[1] > 0.0:
                        if len(sequence) > len_query * length_discrepancy[1]:
                            continue

                    final_dict[query_id].add((name, description, sequence))

    #print("centroid_dict")
    #print(centroid_dict)

    #print("cluster_dict")
    #print(cluster_dict)

    #print("final_dict")
    #print(final_dict)

    # write the fasta files
    print(f"Printing output files...")
    for query_id, hit_list in final_dict.items():
        # update query_id to ensure can generate file name
        query_id_updated = query_id.replace("/", "_")
        query_id_updated = query_id.replace(":", "_")
        query_id_updated = query_id.replace("\\", "_")
        with open(outpref + "_" + query_id + ".fasta", "w") as o:
            for name, description, sequence in hit_list:
                o.write(">" + description + "\n" + sequence + "\n")
    
    sys.exit(0)

                
















if __name__ == "__main__":
    main()