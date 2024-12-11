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
        dir1 = directory(f"{config['output_dir']}/pyrodigal_batches")
        dir2 = directory(f"{config['output_dir']}/pyrodigal_annotations")
        dir3 = directory(f"{config['output_dir']}/mmseqs2_batches")
    shell:
        """
        mkdir {output.dir1}
        mkdir {output.dir2}
        mkdir {output.dir3}
        """

rule split_batch:
    input:
        file_list = f"{config['file_list']}",
        pyrodigal_batch_size = f"{config['pyrodigal_batch_size']}"
    output:
        batch_file = f"{config['output_dir']}/pyrodigal_batches/pyrodigal_batch_N_"
    shell:
        """
        split -l {input.pyrodigal_batch_size} {input.file_list} {output.batch_file}
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
        mmseqs2_batch_size: f"{config['mmseqs2_batch_size']}"
    output:
        pyrodigal_faa_paths: f"{config['output_dir']}/pyrodigal_faa_paths.txt"
        pyrodigal_batches_dir = directory(f"{config['output_dir']}/pyrodigal_batches")
    shell:
        """
        ls -d -1 {input.pyrodigal_annotations}/*/*.faa > {output.pyrodigal_faa_paths}
        split -l {input.mmseqs2_batch_size} {output.pyrodigal_faa_paths} {output.pyrodigal_batches_dir}/pyrodigal_faa_batches_
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
        """


#rule mmseqs2:

