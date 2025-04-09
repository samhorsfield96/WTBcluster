from Bio import SeqIO
import argparse
import pickle

def get_options():
    description = "Annotates tokenised gene clusters"
    parser = argparse.ArgumentParser(description=description,
                                        prog='python annotate_tokens.py')
    IO = parser.add_argument_group('Input/options.out')
    IO.add_argument('--clusters',
                    required=True,
                    help='Cluster file.')
    IO.add_argument('--annotations',
            required=True,
            help='Path to gene annotation tsv file.')
    IO.add_argument('--outpref',
                default="annotated_tokens",
                help='Output prefix.')
    return parser.parse_args()

def main():
    options = get_options()
    clusters = options.clusters
    annotations = options.annotations
    outpref = options.outpref
    
    annotation_dict = {}
    # read in annotations
    with open(annotations, "r") as f:
        # ignore headers
        f.readline()
        f.readline()

        # get column names
        header = f.readline().rstrip()

        #print(header)
        for line in f:
            split_line = line.rstrip().split("\t")
            name = split_line[0]        
            annotation_dict[name] = split_line

    # read in cluster file
    with open(clusters, 'rb') as handle:
        reps_dict, rep_to_token = pickle.load(handle)
    
    # iterate through reps_dict, assigning details to each entry
    with open(outpref + ".tsv", "w") as o:
        o.write("Name\tToken\tAlternate_" + header + "\n")
        for token, name in reps_dict.items():
            annotation = annotation_dict[name]
            annotation_joined = "\t".join(annotation)
            #print(f"{token}, {name}")
            #print(f"{annotation_joined}")

            o.write("{}\t{}\t{}\n".format(str(name), str(token), annotation_joined))


if __name__ == "__main__":
    main()