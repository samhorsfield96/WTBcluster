import glob
import os

configfile: "config.yaml"

rule mkdir:
    output:
        dir1 = directory(f"{config['output_dir']}/pyrodigal_batches"),
        dir2 = directory(f"{config['output_dir']}/pyrodigal_annotations"),
        dir3 = directory(f"{config['output_dir']}/mmseqs2_batches"),
        dir4 = directory(f"{config['output_dir']}/mmseqs2_clustering"),
        dir5 = directory(f"{config['output_dir']}/tokenised_genomes")
    shell:
        """
        mkdir {output.dir1}
        mkdir {output.dir2}
        mkdir {output.dir3}
        mkdir {output.dir4}
        mkdir {output.dir5}
        """

rule split_batch:
    input:
        file_list = f"{config['file_list']}"
    output:
        batch_file = f"{config['output_dir']}/pyrodigal_batches/pyrodigal_batch_N_"
    params:
        pyrodigal_batch_size = f"{config['pyrodigal_batch_size']}"
    shell:
        """
        split -d -l {params.pyrodigal_batch_size} {input.file_list} {output.batch_file}
        """

rule pyrodigal:
    input:
        batch_file = f"{config['output_dir']}/pyrodigal_batches/pyrodigal_batch_N_{{sample}}"
    output:
        outdir = directory(f"{config['output_dir']}/pyrodigal_annotations/batch_{{sample}}_ann")
    conda:
        "WTBcluster"
    threads: 1
    script: "scripts/run_pyrodigal.py"

rule list_pyrodigal:
    input:
        pyrodigal_annotations= f"{config['output_dir']}/pyrodigal_annotations"
    output:
        pyrodigal_faa_paths= f"{config['output_dir']}/pyrodigal_faa_paths.txt",
        pyrodigal_gff_paths= f"{config['output_dir']}/pyrodigal_gff_paths.txt",
        pyrodigal_batches_dir = directory(f"{config['output_dir']}/pyrodigal_batches")
    params:
        mmseqs2_batch_size= f"{config['mmseqs2_batch_size']}"
    shell:
        """
        ls -d -1 {input.pyrodigal_annotations}/*/*.faa > {output.pyrodigal_faa_paths}
        ls -d -1 {input.pyrodigal_annotations}/*/*.gff > {output.pyrodigal_gff_paths}
        split -d -l {params.mmseqs2_batch_size} {output.pyrodigal_faa_paths} {output.pyrodigal_batches_dir}/pyrodigal_faa_batches_
        split -d -l {params.mmseqs2_batch_size} {output.pyrodigal_gff_paths} {output.pyrodigal_batches_dir}/pyrodigal_gff_batches_
        """

rule concatenate_faa:
    input:
        batch_file = f"{config['output_dir']}/pyrodigal_batches/pyrodigal_faa_batches_{{sample}}"
    output:
        outfile = directory(f"{config['output_dir']}/mmseqs2_batches/mmseqs2_concat_batch_{{sample}}.faa")
    conda:
        "WTBcluster"
    threads: 1
    shell:
        """
        > {output.outfile}
        while IFS= read -r line || [[ -n "$line" ]]; do

            cat $line >> {output.outfile}

        done < {input.batch_file}
        """

# Function to list files
def list_files(dir, extension, sort=False):
    input_files = glob.glob(os.path.join(dir, "*" + extension))
    if sort:
        input_files.sort()
    return input_files

# Get number of iterations based on input files
input_files = list_files(f"{config['output_dir']}/mmseqs2_batches", ".faa", sort=True)
num_iterations = len(input_files)

# Define the final output
rule all:
    input:
        expand(f"{config['output_dir']}/mmseqs2_clustering/clustered_{{iteration}}", iteration=range(num_iterations))

# Get current and previous files dynamically
def get_iteration_faa(wildcards):
    return input_files[int(wildcards.iteration)]

def get_iteration_rep(wildcards):
    previous_file = f"{config['output_dir']}/mmseqs2_clustering/clustered_{int(wildcards.iteration) - 1}_rep_seq.fasta"
    if int(wildcards.iteration) == 0 or not os.path.exists(previous_file):
        return None  # No previous file for first iteration
    return previous_file

# Concatenate FASTA files
rule concatenate_fasta:
    input:
        current=lambda wildcards: get_iteration_faa(wildcards),
        previous=lambda wildcards: get_iteration_rep(wildcards)
    output:
        outfile=f"{config['output_dir']}/mmseqs2_clustering/concatenated_{{iteration}}.fasta"
    shell:
        """
        if [ -z "{input.previous}" ]; then
            cp {input.current} {output.outfile}
        else
            cat {input.current} {input.previous} > {output.outfile}
        fi
        """

# Run MMseqs clustering
rule mmseqs_cluster:
    input:
        fasta=f"{config['output_dir']}/mmseqs2_clustering/concatenated_{{iteration}}.fasta"
    output:
        outpref=f"{config['output_dir']}/mmseqs2_clustering/clustered_{{iteration}}"
    threads: 40
    params:
        mmseqs2_tmp_dir=f"{config['mmseqs2_tmp_dir']}",
        mmseqs2_min_ID=f"{config['mmseqs2_min_ID']}",
        mmseqs2_min_cov=f"{config['mmseqs2_min_cov']}",
        mmseqs2_cov_mode=f"{config['mmseqs2_cov_mode']}",
        mmseqs2_ID_mode=f"{config['mmseqs2_ID_mode']}"
    shell:
        """
        mmseqs easy-linclust {input.fasta} {output.outpref} {params.mmseqs2_tmp_dir} \
            --min-seq-id {params.mmseqs2_min_ID} -c {params.mmseqs2_min_cov} \
            --seq-id-mode {params.mmseqs2_ID_mode} --threads {threads} --cov-mode {params.mmseqs2_cov_mode}
        """

rule mmseqs_cluster_merge:
    input:
        indir = directory(f"{config['output_dir']}/mmseqs2_clustering")
    output:
        clusters = f"{config['output_dir']}/mmseqs2_clustering/all_clusters.pkl"
    conda:
        "WTBcluster"
    threads: 1
    script: "scripts/merge_mmseqs2_clusters.py"

rule mmseqs_cluster_write:
    input:
        indir = directory(f"{config['output_dir']}/mmseqs2_clustering"),
        clusters = f"{config['output_dir']}/mmseqs2_clustering/all_clusters.pkl"
    output:
        outfile = f"{config['output_dir']}/mmseqs2_clustering/all_clusters.tsv"
    conda:
        "WTBcluster"
    threads: 1
    script: "scripts/write_mmseqs2_clusters.py"

rule generate_token_db:
    input:
        clusters = f"{config['output_dir']}/mmseqs2_clustering/all_clusters.tsv"
    output:
        out_db = f"{config['output_dir']}/mmseqs2_clustering/gene_tokens.db",
        out_reps = f"{config['output_dir']}/mmseqs2_clustering/reps.pkl"
    conda:
        "WTBcluster"
    threads: 1
    script: "scripts/generate_token_db.py"

rule tokenise_genomes:
    input:
        batch_file = f"{config['output_dir']}/pyrodigal_batches/pyrodigal_gff_batches_{{sample}}",
        out_db = f"{config['output_dir']}/mmseqs2_clustering/gene_tokens.db"
    output:
        outfile= f"{config['output_dir']}/tokenised_genomes/tokenized_genomes_batch_{{sample}}.txt"
    conda:
        "WTBcluster"
    threads: 1
    script: "scripts/tokenize_genomes.py"

# add bakta annotation