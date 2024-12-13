import glob
import os

configfile: "config.yaml"

# Define the final output
rule all:
    input:
        pyrodigal_faa_paths=f"{config['output_dir']}/pyrodigal_faa_paths.txt",
        clusters_pkl = f"{config['output_dir']}/mmseqs2_clustering/all_clusters.pkl",
        clusters_tsv = f"{config['output_dir']}/mmseqs2_clustering/all_clusters.tsv",
        out_db = f"{config['output_dir']}/mmseqs2_clustering/gene_tokens.db",
        out_reps = f"{config['output_dir']}/mmseqs2_clustering/reps.pkl"

# rule mkdir:
#     output:
#         dir1 = directory(f"{config['output_dir']}/pyrodigal_batches"),
#         dir2 = directory(f"{config['output_dir']}/pyrodigal_annotations"),
#         dir3 = directory(f"{config['output_dir']}/mmseqs2_batches"),
#         dir4 = directory(f"{config['output_dir']}/mmseqs2_clustering"),
#         dir5 = directory(f"{config['output_dir']}/tokenised_genomes")
#     shell:
#         """
#         mkdir {output.dir1}
#         mkdir {output.dir2}
#         mkdir {output.dir3}
#         mkdir {output.dir4}
#         mkdir {output.dir5}
#         """

checkpoint split_file_batch:
    input:
        file_list = f"{config['file_list']}"
    output:
        #pyrodigal_batches_dir=directory(f"{config['output_dir']}/pyrodigal_batches"),
        batches_dir=directory(f"{config['output_dir']}/pyrodigal_batches")
    params:
        pyrodigal_batch_size = f"{config['pyrodigal_batch_size']}"
    shell:
        """
        mkdir {output.batches_dir}
        mkdir {output.annotations_dir}
        split -d -l {params.pyrodigal_batch_size} {input.file_list} {output.batches_dir}/pyrodigal_batch_N_
        """

def get_file_batches(wildcards):
    checkpoint_output = checkpoints.split_file_batch.get(wildcards)
    batches_dir = checkpoint_output.output.batches_dir
    return sorted(glob.glob(f"{batches_dir}/pyrodigal_batch_N_*"))

rule pyrodigal:
    input:
        batch_file=lambda wildcards: get_file_batches(wildcards)
    output:
        outdir=directory(f"{config['output_dir']}/pyrodigal_annotations/batch_{{wildcards}}_ann")
    conda:
        "WTBcluster"
    threads: 1
    script: "scripts/run_pyrodigal.py"

checkpoint list_pyrodigal_batch:
    input:
        pyrodigal_annotations= f"{config['output_dir']}/pyrodigal_annotations/batch_{{sample}}_ann"
    output:
        pyrodigal_faa_paths= f"{config['output_dir']}/pyrodigal_batches/pyrodigal_faa_paths_batch_{{sample}}.txt",
        pyrodigal_gff_paths= f"{config['output_dir']}/pyrodigal_batches/pyrodigal_gff_paths_batch_{{sample}}.txt",
    shell:
        """
        ls -d -1 {input.pyrodigal_annotations}/*.faa > {output.pyrodigal_faa_paths}
        ls -d -1 {input.pyrodigal_annotations}/*.gff > {output.pyrodigal_gff_paths}
        """

def get_pyrodigal_batches(wildcards, faa=True):
    checkpoint_output = checkpoints.list_pyrodigal_batch.get(wildcards)
    batches_dir = f"{config['output_dir']}/pyrodigal_batches"
    if faa == True:
        return sorted(glob.glob(f"{batches_dir}/pyrodigal_faa_paths_batch_*.txt"))
    else:
        return sorted(glob.glob(f"{batches_dir}/pyrodigal_gff_paths_batch_*.txt"))

checkpoint list_pyrodigal_full:
    input:
        pyrodigal_faa= lambda wildcards: get_pyrodigal_batches(wildcards, faa=True),
        pyrodigal_gff= lambda wildcards: get_pyrodigal_batches(wildcards, faa=False)
    output:
        pyrodigal_faa_paths= f"{config['output_dir']}/pyrodigal_batches/pyrodigal_faa_paths_all.txt",
        pyrodigal_gff_paths= f"{config['output_dir']}/pyrodigal_batches/pyrodigal_gff_paths_all.txt",
        pyrodigal_batches_dir= f"{config['output_dir']}/pyrodigal_batches",
    params:
        mmseqs2_batch_size= f"{config['mmseqs2_batch_size']}"
    shell:
        """
        cat {input.pyrodigal_faa} > {output.pyrodigal_faa_paths}
        cat {input.pyrodigal_gff} > {output.pyrodigal_gff_paths}
        split -d -l {params.mmseqs2_batch_size} {output.pyrodigal_faa_paths} {output.pyrodigal_batches_dir}/pyrodigal_faa_batches_
        split -d -l {params.mmseqs2_batch_size} {output.pyrodigal_gff_paths} {output.pyrodigal_batches_dir}/pyrodigal_gff_batches_
        """

def get_mmseqs2_batches(wildcards):
    checkpoint_output = checkpoints.list_pyrodigal_full.get(wildcards)
    batches_dir = checkpoint_output.output.pyrodigal_batches_dir
    return sorted(glob.glob(f"{batches_dir}/pyrodigal_faa_batches_*"))

rule concatenate_faa:
    input:
        batch_file=lambda wildcards: get_mmseqs2_batches(wildcards)
    output:
        outfile = directory(f"{config['output_dir']}/mmseqs2_batches/mmseqs2_concat_batch_{{wildcards}}.faa")
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

rule iteration:
    input:
        expand(f"{config['output_dir']}/mmseqs2_clustering/clustered_{{iteration}}", iteration=range(num_iterations)),

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
        indir = f"{config['output_dir']}/mmseqs2_clustering"
    output:
        clusters = f"{config['output_dir']}/mmseqs2_clustering/all_clusters.pkl"
    conda:
        "WTBcluster"
    threads: 1
    script: "scripts/merge_mmseqs2_clusters.py"

rule mmseqs_cluster_write:
    input:
        indir = f"{config['output_dir']}/mmseqs2_clustering",
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