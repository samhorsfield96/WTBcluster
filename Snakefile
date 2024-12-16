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
        f"{config['output_dir']}/list_pyrodigal_full.done"
        #f"{config['output_dir']}/concatenate_faa.done",
        #f"{config['output_dir']}/tokenisation.done",
        #out_db = f"{config['output_dir']}/mmseqs2_clustering/gene_tokens.db",

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

checkpoint list_pyrodigal_full:
    input:
        dir_list=get_pyrodigal_outputs,
    output:
        check_file=f"{config['output_dir']}/list_pyrodigal_full.done"
    params:
        mmseqs2_num_batches=f"{config['mmseqs2_num_batches']}",
        batches_dir=f"{config['output_dir']}/pyrodigal_batches",
        outpref="mmseqs2_batch_N"
    script: "scripts/split_file_mmseqs.py"

# def get_mmseqs2_batches(wildcards):
#     checkpoint_output = checkpoints.list_pyrodigal_full.get(**wildcards).output[0]
#     return expand("{out_dir}/pyrodigal_batches/mmseqs2_batch_N{batch_ID}_faa.txt", out_dir=config['output_dir'], batch_ID=range(config['mmseqs2_num_batches']))

# rule check_list_pyrodigal_full:
#     input:
#         get_mmseqs2_batches
#     output:
#         touch(f"{config['output_dir']}/list_pyrodigal_full.done")
#     run:
#         pass

# rule concatenate_faa:
#     input:
#         batch_file = f"{config['output_dir']}/pyrodigal_batches/mmseqs2_batch_N{{batch_ID}}_faa.txt",
#     output:
#         outfile = f"{config['output_dir']}/mmseqs2_batches/mmseqs2_concat_batch_N{{batch_ID}}.faa"
#     conda:
#         "WTBcluster"
#     threads: 1
#     shell:
#         """
#         > {output.outfile}
#         while IFS= read -r line || [[ -n "$line" ]]; do

#             cat $line >> {output.outfile}

#         done < {input.batch_file}
#         """

# def get_mmseqs2_concat_batches(wildcards):
#     checkpoint_output = checkpoints.concatenate_faa.get(**wildcards).output[0]
#     return expand("{out_dir}/mmseqs2_batches/mmseqs2_concat_batch_N{batch_ID}.faa", out_dir=config['output_dir'], batch_ID=range(config['mmseqs2_num_batches']))

# rule check_concatenate_faa:
#     input:
#         get_mmseqs2_concat_batches
#     output:
#         touch(f"{config['output_dir']}/concatenate_faa.done")
#     run:
#         pass

# # Function to list files
# def list_files(dir, extension, sort=False):
#     input_files = glob.glob(os.path.join(dir, "*" + extension))
#     if sort:
#         input_files.sort()
#     return input_files

# # Get number of iterations based on input files
# input_files = list_files(f"{config['output_dir']}/mmseqs2_batches", ".faa", sort=True)
# num_iterations = len(input_files)

# rule iteration1:
#     input:
#         expand( f"{config['output_dir']}/mmseqs2_batches/mmseqs2_concat_batch_{{iteration}}.faa", iteration=range(num_iterations)),

# # Get current and previous files dynamically
# def get_iteration(wildcards):
#     return input_files[int(wildcards.iteration)]

# def get_iteration_rep(wildcards):
#     previous_file = f"{config['output_dir']}/mmseqs2_clustering/clustered_{int(wildcards.iteration) - 1}_rep_seq.fasta"
#     if int(wildcards.batch_ID) == 0 or not os.path.exists(previous_file):
#         return None  # No previous file for first iteration
#     return previous_file

# # Concatenate FASTA files
# rule concatenate_fasta:
#     input:
#         current=lambda wildcards: get_iteration(wildcards),
#         previous=lambda wildcards: get_iteration_rep(wildcards)
#     output:
#         outfile=f"{config['output_dir']}/mmseqs2_clustering/concatenated_{{iteration}}.fasta"
#     shell:
#         """
#         if [ -z "{input.previous}" ]; then
#             cp {input.current} {output.outfile}
#         else
#             cat {input.current} {input.previous} > {output.outfile}
#         fi
#         """

# # Run MMseqs clustering
# rule mmseqs_cluster:
#     input:
#         fasta=f"{config['output_dir']}/mmseqs2_clustering/concatenated_{{iteration}}.fasta"
#     output:
#         outpref=f"{config['output_dir']}/mmseqs2_clustering/clustered_{{iteration}}"
#     threads: 40
#     params:
#         mmseqs2_tmp_dir=f"{config['mmseqs2_tmp_dir']}",
#         mmseqs2_min_ID=f"{config['mmseqs2_min_ID']}",
#         mmseqs2_min_cov=f"{config['mmseqs2_min_cov']}",
#         mmseqs2_cov_mode=f"{config['mmseqs2_cov_mode']}",
#         mmseqs2_ID_mode=f"{config['mmseqs2_ID_mode']}"
#     shell:
#         """
#         mmseqs easy-linclust {input.fasta} {output.outpref} {params.mmseqs2_tmp_dir} \
#             --min-seq-id {params.mmseqs2_min_ID} -c {params.mmseqs2_min_cov} \
#             --seq-id-mode {params.mmseqs2_ID_mode} --threads {threads} --cov-mode {params.mmseqs2_cov_mode}
#         """

# def get_mmseqs2_clusters(indir):
#     list_of_samples = [path.split("/clustered_")[-1].replace("_cluster.tsv", "") for path in glob.glob(indir + "/clustered_*_cluster.tsv")]
#     return sorted(list_of_samples)

# rule mmseqs_cluster_merge:
#     input:
#         infiles= expand(f"{config['output_dir']}/mmseqs2_clustering/clustered_{{batch_ID}}_cluster.tsv", batch_ID=get_mmseqs2_clusters(f"{config['output_dir']}/mmseqs2_clustering"))
#     output:
#         clusters = f"{config['output_dir']}/mmseqs2_clustering/all_clusters.pkl"
#     conda:
#         "WTBcluster"
#     threads: 1
#     script: "scripts/merge_mmseqs2_clusters.py"

# rule mmseqs_cluster_write:
#     input:
#         infiles = lambda wildcards: get_mmseqs2_clusters(os.path.join(config['output_dir'], "mmseqs2_clustering")),
#         clusters = f"{config['output_dir']}/mmseqs2_clustering/all_clusters.pkl"
#     output:
#         outfile = f"{config['output_dir']}/mmseqs2_clustering/all_clusters.tsv"
#     conda:
#         "WTBcluster"
#     threads: 1
#     script: "scripts/write_mmseqs2_clusters.py"

# rule generate_token_db:
#     input:
#         clusters = f"{config['output_dir']}/mmseqs2_clustering/all_clusters.tsv"
#     output:
#         out_db = f"{config['output_dir']}/mmseqs2_clustering/gene_tokens.db",
#         out_reps = f"{config['output_dir']}/mmseqs2_clustering/reps.pkl"
#     conda:
#         "WTBcluster"
#     threads: 1
#     script: "scripts/generate_token_db.py"

# rule tokenise_genomes:
#     input:
#         batch_file = f"{config['output_dir']}/pyrodigal_batches/pyrodigal_gff_paths_batch_{{batch_ID}}.txt",
#         out_db = f"{config['output_dir']}/mmseqs2_clustering/gene_tokens.db"
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