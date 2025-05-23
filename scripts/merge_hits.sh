#!/bin/bash

# Usage message
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 INPUT_DIR [OUTPUT_DIR]"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="${2:-concatenated_fastas}"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Move to input directory
cd "$INPUT_DIR" || exit 1

# Get unique suffixes after 'batch_X_'
for suffix in $(ls batch_*_*.fasta 2>/dev/null | sed -E 's/batch_[^_]+_//' | sort -u); do
    echo "Concatenating files for $suffix"
    cat batch_*_"$suffix" > "$OUTPUT_DIR/$suffix"
done

echo "Done. Concatenated files are in '$OUTPUT_DIR'."