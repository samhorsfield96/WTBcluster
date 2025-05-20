#!/bin/bash

input_list=/path/to/faa_file_lists.txt
outdir=/path/to/outdir
num_batches=10


# Split into 10 line-balanced batches with numeric suffixes
split -n l/${num_batches} -d "$input_list" ${outdir}/batch_faa_ --additional-suffix=.txt

# Iterate over all batch_* files in sorted order
for batch_file in $(ls ${outdir}/batch_faa_* | sort); do
  # Extract just the numeric suffix, remove leading zeros
    batch_num=$(basename "$batch_file" .txt | sed -E 's/^batch_faa_//' | sed 's/^0*//')
    batch_num=${batch_num:-0}  # fallback in case we get empty string (e.g., "00")

    output_file="${outdir}/batch_${batch_num}.fasta"
    cat $(<"$batch_file") > "$output_file"
done

# Optional: cleanup intermediate files
rm ${outdir}/batch_faa_*.txt