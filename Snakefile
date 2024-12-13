import glob
import os

configfile: "config.yaml"

# Define the final output
rule all:
   input:
        #pyrodigal_batch_dir = f"{config['output_dir']}/pyrodigal_batches",
        #pyrodigal_batch_file = expand(f"{config['output_dir']}/pyrodigal_batches/pyrodigal_batch_N{{batch_ID}}", batch_ID=get_input_batches(f"{config['output_dir']}/pyrodigal_batches")),
        pyrodigal_faa_paths = f"{config['output_dir']}/pyrodigal_batches/pyrodigal_faa_paths_all.txt",
        #pyrodigal_gff_paths = f"{config['output_dir']}/pyrodigal_batches/pyrodigal_gff_paths_all.txt",
        #clusters_pkl = f"{config['output_dir']}/mmseqs2_clustering/all_clusters.pkl",
        #clusters_tsv = f"{config['output_dir']}/mmseqs2_clustering/all_clusters.tsv",
        #out_db = f"{config['output_dir']}/mmseqs2_clustering/gene_tokens.db",
        #out_reps = f"{config['output_dir']}/mmseqs2_clustering/reps.pkl"

rule split_file_batch:
    input:
        file_list=f"{config['file_list']}"
    output:
        batches_dir=directory(f"{config['output_dir']}/pyrodigal_batches")
    params:
        pyrodigal_batch_size=int(f"{config['pyrodigal_batch_size']}"),
        outpref="pyrodigal_batch_N"
    script: "scripts/split_file.py"

rule pyrodigal:
    input:
        batch_file=f"{config['output_dir']}/pyrodigal_batches/pyrodigal_batch_N{{batch_ID}}.txt"
    output:
        outdir=directory(f"{config['output_dir']}/pyrodigal_annotations/batch_{{batch_ID}}_ann")
    conda:
        "WTBcluster"
    threads: 1
    script: "scripts/run_pyrodigal.py"

def get_pyrodigal_batchs(indir):
    list_of_samples = [path.split("/batch_")[-1].replace("_ann", "") for path in glob.glob(indir + "/batch_*_ann")]
    return sorted(list_of_samples)

rule list_pyrodigal_full:
    input:
        pyrodigal_annotations=f"{config['output_dir']}/pyrodigal_annotations/batch_{{batch_ID}}_ann",
    output:
        pyrodigal_faa_paths=f"{config['output_dir']}/pyrodigal_batches/pyrodigal_faa_paths_all.txt",
        pyrodigal_gff_paths=f"{config['output_dir']}/pyrodigal_batches/pyrodigal_gff_paths_all.txt",
    params:
        mmseqs2_batch_size=f"{config['mmseqs2_batch_size']}",
        pyrodigal_batches_dir=f"{config['output_dir']}/pyrodigal_batches"
    shell:
        """
        cat {input.pyrodigal_annotations}/*.faa > {output.pyrodigal_faa_paths}
        cat {input.pyrodigal_annotations}/*.ffn > {output.pyrodigal_gff_paths}
        split -d -l {params.mmseqs2_batch_size} {output.pyrodigal_faa_paths} {params.pyrodigal_batches_dir}/pyrodigal_faa_batches_N
        split -d -l {params.mmseqs2_batch_size} {output.pyrodigal_gff_paths} {params.pyrodigal_batches_dir}/pyrodigal_gff_batches_N
        """

# def get_mmseqs2_batches(indir, faa=True):
#     list_of_samples = [path.split("_N")[-1] for path in glob.glob(indir + "/pyrodigal_faa_batches_N*")]
#     return sorted(list_of_samples)

# rule concatenate_faa:
#     input:
#         batch_file = expand(f"{config['output_dir']}/pyrodigal_batches/pyrodigal_faa_batches_N{{batch_ID}}", batch_ID=get_mmseqs2_batches(f"{config['output_dir']}/pyrodigal_batches")),
#     output:
#         outfile = f"{config['output_dir']}/mmseqs2_batches/mmseqs2_concat_batch_{{batch_ID}}.faa"
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
#         infiles = lambda wildcards: get_mmseqs2_clusters(wildcards),
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

# add bakta annotation