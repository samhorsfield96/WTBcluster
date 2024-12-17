import glob
import os

configfile: "config.yaml"

def get_pyrodigal_outputs(wildcards):
    checkpoint_output = checkpoints.split_file_batch.get(**wildcards).output[0]
    return expand("{out_dir}/pyrodigal_annotations/batch_{batch_ID}_ann", out_dir=config['output_dir'], batch_ID=range(config['pyrodigal_num_batches']))

def get_tokenisation_outputs(wildcards):
    return expand("{out_dir}/tokenised_genomes/tokenized_genomes_batch_{batch_ID}.txt", out_dir=config['output_dir'], batch_ID=range(config['mmseqs2_num_batches']))
    
# Define the final output
rule all:
   input:
        #f"{config['output_dir']}/pyrodigal.done",
        f"{config['output_dir']}/list_pyrodigal_full.done",
        #f"{config['output_dir']}/concatenate_faa.done",
        #f"{config['output_dir']}/mmseqs_cluster.done",
        f"{config['output_dir']}/tokenised_genomes/all_clusters.pkl",
        f"{config['output_dir']}/tokenised_genomes/all_clusters.tsv",
        f"{config['output_dir']}/tokenised_genomes/reps.pkl",
        #f"{config['output_dir']}/tokenisation.done",

checkpoint split_file_batch:
    input:
        file_list=f"{config['file_list']}",
    output:
        batches_dir=directory(f"{config['output_dir']}/pyrodigal_batches")
    params:
        pyrodigal_num_batches=int(f"{config['pyrodigal_num_batches']}"),
        outpref="pyrodigal_batch_N"
    script: "scripts/split_file_pyrodigal.py"

checkpoint pyrodigal:
    input:
        batch_file=f"{config['output_dir']}/pyrodigal_batches/pyrodigal_batch_N{{batch_ID}}.txt"
    output:
        outdir=directory(f"{config['output_dir']}/pyrodigal_annotations/batch_{{batch_ID}}_ann"),
    conda:
        "WTBcluster"
    threads: 1
    script: "scripts/run_pyrodigal.py"

# rule check_pyrodigal:
#     input:
#         get_pyrodigal_outputs
#     output:
#         touch(f"{config['output_dir']}/pyrodigal.done")
#     run:
#         pass

rule list_pyrodigal_full:
    input:
        dir_list=get_pyrodigal_outputs,
    output:
        check_file=f"{config['output_dir']}/list_pyrodigal_full.done"
    params:
        mmseqs2_num_batches=f"{config['mmseqs2_num_batches']}",
        batches_dir=f"{config['output_dir']}/pyrodigal_batches",
        outpref="mmseqs2_batch_N"
    script: "scripts/split_file_mmseqs.py"

checkpoint concatenate_faa:
    input:
        batch_file = f"{config['output_dir']}/pyrodigal_batches/mmseqs2_batch_N{{batch_ID}}_faa.txt",
    output:
        outfile = f"{config['output_dir']}/mmseqs2_batches/mmseqs2_concat_batch_N{{batch_ID}}.faa"
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

def get_concatentated_files():
    return expand("{out_dir}/mmseqs2_batches/mmseqs2_concat_batch_N{batch_ID}.faa", out_dir=config['output_dir'], batch_ID=range(config['mmseqs2_num_batches']))

# Run MMseqs clustering
checkpoint mmseqs_cluster:
    input:
        file_list=get_concatentated_files()
    output:
        output_dir=directory(f"{config['output_dir']}/mmseqs2_clustering")
    threads: 40
    params:
        mmseqs2_tmp_dir=f"{config['mmseqs2_tmp_dir']}",
        mmseqs2_min_ID=f"{config['mmseqs2_min_ID']}",
        mmseqs2_min_cov=f"{config['mmseqs2_min_cov']}",
        mmseqs2_cov_mode=f"{config['mmseqs2_cov_mode']}",
        mmseqs2_ID_mode=f"{config['mmseqs2_ID_mode']}",
        outpref="clustering_",
        
    script: "scripts/run_mmseqs2.py"

def get_mmseqs2_clusters(wildcards):
    checkpoint_output = checkpoints.mmseqs_cluster.get(**wildcards).output[0]
    return expand("{out_dir}/mmseqs2_clustering/clustering_{batch_ID}_cluster.tsv", out_dir=config['output_dir'], batch_ID=range(config['mmseqs2_num_batches']))

# rule check_mmseqs_cluster:
#     input:
#         get_mmseqs2_clusters
#     output:
#         touch(f"{config['output_dir']}/mmseqs_cluster.done")
#     run:
#         pass

rule mmseqs_cluster_merge:
    input:
        infiles=get_mmseqs2_clusters
    output:
        clusters = f"{config['output_dir']}/tokenised_genomes/all_clusters.pkl"
    conda:
        "WTBcluster"
    threads: 1
    script: "scripts/merge_mmseqs2_clusters.py"

rule mmseqs_cluster_write:
    input:
        infiles=get_mmseqs2_clusters,
        clusters = f"{config['output_dir']}/tokenised_genomes/all_clusters.pkl"
    output:
        outfile = f"{config['output_dir']}/tokenised_genomes/all_clusters.tsv"
    conda:
        "WTBcluster"
    threads: 1
    script: "scripts/write_mmseqs2_clusters.py"

rule generate_token_db:
    input:
        clusters = f"{config['output_dir']}/tokenised_genomes/all_clusters.tsv"
    output:
        out_db = directory(f"{config['output_dir']}/tokenised_genomes/gene_tokens.db"),
        out_reps = f"{config['output_dir']}/tokenised_genomes/reps.pkl"
    conda:
        "WTBcluster"
    threads: 1
    script: "scripts/generate_token_db.py"

# rule tokenise_genomes:
#     input:
#         batch_file = f"{config['output_dir']}/pyrodigal_batches/pyrodigal_gff_paths_batch_{{batch_ID}}.txt",
#         out_db = f"{config['output_dir']}/tokenised_genomes/gene_tokens.db"
#     output:
#         outfile= f"{config['output_dir']}/tokenised_genomes/tokenized_genomes_batch_{{batch_ID}}.txt",
#     conda:
#         "WTBcluster"
#     threads: 1
#     script: "scripts/tokenize_genomes.py"

# rule aggregate_tokens:
#     input:
#         get_tokenisation_outputs
#     output:
#         touch(f"{config['output_dir']}/tokenisation.done")
#     run:
#         pass

# add bakta annotation