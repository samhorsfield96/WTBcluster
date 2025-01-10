from Bio import SeqIO
import pyrodigal
import gzip
import os
import subprocess
from natsort import natsorted

def run_mmseqs(file_list, output_dir, outpref, mmseqs2_tmp_dir, mmseqs2_min_ID, mmseqs2_min_cov, mmseqs2_cov_mode, mmseqs2_cluster_mode, mmseqs2_ID_mode, threads):
    """Trains and runs mmseqs2 iteratively."""
    
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # sort concatenated input files
    file_list = natsorted(file_list)
    
    # determine iteration number
    current_index = 0
    for index, filename in enumerate(file_list):
        rep_file_path = f"{output_dir}/{outpref}{index}_rep_seq.fasta"
        if not os.path.exists(rep_file_path):
            current_index = index
            break
    
    # check to see if all files complete
    if current_index < len(file_list) - 1:
        # start by getting startign file
        for index in range(current_index, len(file_list)):
            current_file = file_list[index]

            mmseqs_input_file = os.path.join(output_dir, f"concatenated_{index}.fasta")

            # concatenate representative 
            if index != 0:
                prev_file = f"{output_dir}/{outpref}{index - 1}_rep_seq.fasta"

                # concatenate files
                print(f"Concatentating {current_file} and {prev_file} to {mmseqs_input_file}")
                with open(mmseqs_input_file, "w") as o:
                    with open(current_file, "r") as f1:
                        o.write(f1.read())
                    with open(prev_file, "r") as f2:
                        o.write(f2.read())

            else:
                # copy original file
                print(f"Copying {current_file} to {mmseqs_input_file}")
                with open(mmseqs_input_file, "w") as o:
                    with open(current_file, "r") as f1:
                        o.write(f1.read())

            output_prefix = os.path.join(output_dir, outpref + str(index))

            try:
                # Step 2: Run linclust with parameters
                subprocess.run([
                    "mmseqs", "easy-linclust", mmseqs_input_file, output_prefix, mmseqs2_tmp_dir,
                    "--min-seq-id", str(mmseqs2_min_ID),
                    "-c", str(mmseqs2_min_cov),
                    "--seq-id-mode", str(mmseqs2_ID_mode),
                    "--cov-mode", str(mmseqs2_cov_mode),
                    "--cluster-mode", str(mmseqs2_cluster_mode),
                    "--threads", str(threads),
                    "-v", str(2) # only print errors and warnings
                ], check=True)

                print(f"Clustering results saved to {output_prefix}")

            except subprocess.CalledProcessError as e:
                print(f"An error occurred: {e}")
        
run_mmseqs(snakemake.input.file_list, snakemake.output.output_dir, snakemake.params.outpref, snakemake.params.mmseqs2_tmp_dir, snakemake.params.mmseqs2_min_ID, snakemake.params.mmseqs2_min_cov, snakemake.params.mmseqs2_cov_mode, snakemake.params.mmseqs2_cluster_mode, snakemake.params.mmseqs2_ID_mode, snakemake.threads)