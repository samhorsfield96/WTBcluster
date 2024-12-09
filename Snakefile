import glob

configfile: "config.yaml"

rule all:
    input:
        #matrix = f"{config['output_dir']}/presence_absence_matrix.txt",
        #summary_file = f"{config['output_dir']}/pangenome_summary.tsv",
        #checkm_file = f"{config['output_dir']}/checkm_out.tsv",
        #cgt_output = f"{config['output_dir']}/cgt_output.txt"

rule mkdir:
    output:
        dir1 = directory(f"{config['output_dir']}/prodigal_batches")
        dir2 = directory(f"{config['output_dir']}/prodigal_annotations")
    shell:
        """
        mkdir {output.dir1}
        mkdir {output.dir2}
        """

rule split_batch:
    input:
        file_list = f"{config['file_list']}",
        batch_size = f"{config['batch_size']}"
    output:
        batch_file = f"{config['output_dir']}/prodigal_batches/batch_N_"
    shell:
        """
        split -l {input.batch_size} {input.file_list} {output.batch_file}
        """

rule pyrodigal:
    input:
        batch_file = f"{config['output_dir']}/prodigal_batches/batch_N_{{sample}}"
    output:
        outdir = directory(f"{config['output_dir']}/prodigal_annotations/batch_{{sample}}_ann")
    conda:
        "WTBcluster"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 15000
    script: "scripts/run_pyrodigal.py"
