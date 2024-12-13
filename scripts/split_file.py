import os
import sys

def split_file(input_file, lines_per_file, output_dir, outpref):
    """
    Splits a file into multiple smaller files with a specified number of lines each.

    Parameters:
    input_file (str): Path to the input file.
    lines_per_file (int): Number of lines per output file.
    output_dir (str): Directory to save the output files.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    try:
        with open(input_file, 'r') as infile:
            file_count = 0
            lines_written = 0
            output_file = None

            for line in infile:
                if lines_written % lines_per_file == 0:
                    if output_file:
                        output_file.close()
                    file_count += 1
                    output_path = os.path.join(output_dir, f"{outpref}{file_count}.txt")
                    output_file = open(output_path, 'w')

                output_file.write(line)
                lines_written += 1

            if output_file:
                output_file.close()

        print(f"Split complete: {file_count} files created in '{output_dir}'.")

    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

split_file(snakemake.input.file_list, snakemake.params.pyrodigal_batch_size, snakemake.output.batches_dir, snakemake.params.outpref)