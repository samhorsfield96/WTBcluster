from Bio import SeqIO
import pyrodigal
import gzip
import os
import subprocess
from natsort import natsorted

def run_mmseqs(file_list, output_dir, mmseqs2_tmp_dir, mmseqs2_min_ID, mmseqs2_min_cov, mmseqs2_cov_mode, mmseqs2_cluster_mode, mmseqs2_ID_mode, mmseqs2_alignment_mode, threads):
    """Merges parallel MMseqs2 runs."""
    
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # sort concatenated input files
    file_list = natsorted(file_list)
    
    # create final file for clustering
    mmseqs_input_file = os.path.join(output_dir, f"concatenated_final.fasta")
    with open(mmseqs_input_file, "w") as o:
            for current_file in file_list:
                print(f"Concatentating {current_file} to {mmseqs_input_file}")
                with open(current_file, "r") as f1:
                    while True:
                        line = f1.readline()
                        if not line:
                            break
                        o.write(line)

    output_prefix = os.path.join(output_dir, "final")

    try:
        # Step 2: Run linclust with parameters
        subprocess.run([
            "mmseqs", "easy-linclust", mmseqs_input_file, output_prefix, mmseqs2_tmp_dir,
            "--min-seq-id", str(mmseqs2_min_ID),
            "-c", str(mmseqs2_min_cov),
            "--seq-id-mode", str(mmseqs2_ID_mode),
            "--cov-mode", str(mmseqs2_cov_mode),
            "--cluster-mode", str(mmseqs2_cluster_mode),
            "--alignment-mode" , str(mmseqs2_alignment_mode),
            "--threads", str(threads),
            "-v", str(2) # only print errors and warnings
        ], check=True)

        print(f"Clustering results saved to {output_prefix}")

    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}")
        
run_mmseqs(snakemake.input.file_list, snakemake.output.output_dir, snakemake.params.mmseqs2_tmp_dir, snakemake.params.mmseqs2_min_ID, snakemake.params.mmseqs2_min_cov, snakemake.params.mmseqs2_cov_mode, snakemake.params.mmseqs2_cluster_mode, snakemake.params.mmseqs2_ID_mode, snakemake.params.mmseqs2_alignment_mode, snakemake.threads)