import glob
import os

configfile: "config.yaml"

rule all:
    input:
        #matrix = f"{config['output_dir']}/presence_absence_matrix.txt",
        #summary_file = f"{config['output_dir']}/pangenome_summary.tsv",
        #checkm_file = f"{config['output_dir']}/checkm_out.tsv",
        #cgt_output = f"{config['output_dir']}/cgt_output.txt"

rule mkdir:
    output:
        dir1 = directory(f"{config['output_dir']}/pyrodigal_batches"),
        dir2 = directory(f"{config['output_dir']}/pyrodigal_annotations"),
        dir3 = directory(f"{config['output_dir']}/mmseqs2_batches")
        dir4 = directory(f"{config['output_dir']}/mmseqs2_clustering")
    shell:
        """
        mkdir {output.dir1}
        mkdir {output.dir2}
        mkdir {output.dir3}
        mkdir {output.dir4}
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
        pyrodigal_annotations: f"{config['output_dir']}/pyrodigal_annotations"
    output:
        pyrodigal_faa_paths: f"{config['output_dir']}/pyrodigal_faa_paths.txt",
        pyrodigal_batches_dir = directory(f"{config['output_dir']}/pyrodigal_batches")
    params:
        mmseqs2_batch_size: f"{config['mmseqs2_batch_size']}"
    shell:
        """
        ls -d -1 {input.pyrodigal_annotations}/*/*.faa > {output.pyrodigal_faa_paths}
        split -d -l {params.mmseqs2_batch_size} {output.pyrodigal_faa_paths} {output.pyrodigal_batches_dir}/pyrodigal_faa_batches_
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

# get list of input files
def list_files(dir, extension, sort=False):
    input_files = glob.glob(os.path.join(dir, "*" + extension))
    if sort:
        input_files.sort()
    return input_files

def get_iteration_faa(wildcards, dir, extension):
    input_files = list_files(dir, extension, sort=True)
    iteration = int(wildcards.iteration)
    return input_files[iteration]

def get_iteration_rep(wildcards, dir, extension):
    input_files = list_files(dir, extension, sort=True)
    iteration = int(wildcards.iteration)
    if iteration == 0:
        return None
    return input_files[iteration]


rule concatenate_fasta:
    input:
        current=lambda wildcards: expand(f"{{sample}}", sample=get_iteration_faa(wildcards, f"{config['output_dir']}/mmseqs2_batches", ".faa")),
        previous=lambda wildcards: expand(f"{{sample}}", sample=get_iteration_rep(wildcards, f"{config['output_dir']}/mmseqs2_clustering", "_rep_seq.fasta"))
    output:
        outfile: f"{config['output_dir']}/mmseqs2_clustering/concatenated_{iteration}.fasta"
    shell:
        """
        cat {input.current} {input.previous} > {output.outfile} || cp {input.current} {output}
        """

rule mmseqs_cluster:
    input:
        fasta= f"{config['output_dir']}/mmseqs2_clustering/concatenated_{iteration}.fasta"
    output:
        outpref=f"{config['output_dir']}/mmseqs2_clustering/clustered_{iteration}"
    threads: 40
    params:
        mmseqs2_tmp_dir=f"{config['mmseqs2_tmp_dir']}",
        mmseqs2_min_ID= f"{config['mmseqs2_min_ID']}",
        mmseqs2_min_cov= f"{config['mmseqs2_min_cov']}",
        mmseqs2_cov_mode= f"{config['mmseqs2_cov_mode']}",
        mmseqs2_ID_mode= f"{config['mmseqs2_ID_mode']}",
    shell:
        """
        mmseqs easy-linclust {input.fasta} {output.outpref} {params.mmseqs2_tmp_dir} --min-seq-id {params.mmseqs2_min_ID} -c {params.mmseqs2_min_cov} --seq-id-mode {params.mmseqs2_ID_mode} --threads {threads} --cov-mode {params.mmseqs2_cov_mode}
        """

