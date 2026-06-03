import pandas as pd
import argparse
from scipy.sparse import lil_matrix
import csv
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description="Generate gene pa matrix from tokenised genomes.py.")
    parser.add_argument("--infile", type=str, required=True, help="Path to tokenised genome file.")
    parser.add_argument("--outpref", default="output", help="Output prefix")
    parser.add_argument("--reverse", default=False, action="store_true", help="Reverse output so each row is genome and each column is gene")

    args = parser.parse_args()

    return args

def main():
    options = parse_args()
    infile = options.infile
    outpref = options.outpref
    reverse = options.reverse

    # Read file
    genome_list = []
    genome_index = {}

    # first pass through file, determine number of genes
    print("Counting data size...")
    total_genes_set = set()
    total_genomes = 0
    with open(infile, "r") as f:
        for line in f:
            genome, vals = line.strip().split("\t")
            vals = vals.split()

            vals = [abs(int(v)) for v in vals if v != "_"]

            total_genes_set.update(vals)
            genome_index[genome] = total_genomes
            genome_list.append(genome)
            total_genomes += 1
            
    # total_genes_set is your input set
    sorted_genes = sorted(total_genes_set)

    # Create an index mapping: gene_id → row_index
    gene_index = {gene: i for i, gene in enumerate(sorted_genes)}

    total_genes = len(total_genes_set)
    print(f"total_genes: {total_genes}")
    print(f"total_genomes: {total_genomes}")

    if reverse:
        M = lil_matrix((total_genomes, total_genes), dtype=bool)
    else:
        M = lil_matrix((total_genes, total_genomes), dtype=bool)

    print("Filling pa matrix...")
    with open(infile, "r") as f:
        for line in f:
            genome, vals = line.strip().split("\t")
            vals = vals.split()

            vals = [abs(int(v)) for v in vals if v != "_"]

            genome_ID = genome_index[genome]
            for val in vals:
                gene_ID = gene_index[val]
                if reverse:
                    M[genome_ID, gene_ID] = 1
                else:
                    M[gene_ID, genome_ID] = 1

    print("Writing file...")
    with open(outpref + "_pa.csv", "w", newline="") as f:
        writer = csv.writer(f)

        # Optional: write header
        if reverse:
            writer.writerow(["genome"] + sorted_genes)
        else:
            writer.writerow(["gene"] + genome_list)

        for i in range(M.shape[0]):
            # Get sparse row as coordinate lists
            row = M.rows[i]      # list of column indices where value = True
            vals = M.data[i]     # list of True values (all True for bool matrix)

            # Build dense row efficiently
            dense = np.zeros(M.shape[1], dtype=int)
            dense[row] = 1       # mark present positions

            if reverse:
                writer.writerow([genome_list[i]] + dense.tolist())
            else:
                writer.writerow([sorted_genes[i]] + dense.tolist())

    # write index files
    with open(outpref + "_gene_index.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["gene_id", "index"])
        for gene, idx in gene_index.items():
            writer.writerow([gene, idx])

    with open(outpref + "_genome_index.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["genome_id", "index"])
        for genome, idx in genome_index.items():
            writer.writerow([genome, idx])

if __name__ == "__main__":
    main()
