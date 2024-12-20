# WTBcluster: a snakemake pipeline for clustering loads of proteins

WTBcluster (stands for 'Woooowww That's Big-cluster') calls bacterial proteins using [Prodigal](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-119) iteratively clusters proteins using [Linclust](https://www.nature.com/articles/s41467-018-04964-5), part of the [MMseqs2](https://www.nature.com/articles/nbt.3988) suite of tools.

## Dependencies:

* python>=3.9
* pyrodigal
* mmseqs2
* pandas
* biopython
* snakemake
* python-rocksdb
* natsort

## To install:

Install the required packages using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)/[mamba](https://github.com/mamba-org/mamba):

```
git clone https://github.com/samhorsfield96/WTBcluster.git
cd WTBcluster
mamba env create -f environment.yml
mamba activate WTBcluster
```

# Quick start:

Update `config.yaml` to specify workflow and directory paths.
- `output_dir`: path to output directory. Does not need to exist prior to running.
- `file_list`: path to file containing paths to FASTA files (gzipped or uncompressed), one per line.
- `pyrodigal_num_batches`: number of pyrodigal batches to split the `file_list` into.
- `mmseqs2_num_batches`: number of mmseqs2 batches to split the `file_list` into.

- `mmseqs2_min_ID`: minimum identity between two CDS to be clustered.
- `mmseqs2_min_cov`: minimum sequence coverage between two CDS to be clustered.
- `mmseqs2_cov_mode`: method to calculate coverage between two CDS (See [guide](https://mmseqs.com/latest/userguide.pdf))
- `mmseqs2_ID_mode`: method to calculate identity between two CDS (See [guide](https://mmseqs.com/latest/userguide.pdf))
- `mmseqs2_tmp_dir`: temporary directory for MMseqs2 processing.


Run snakemake (must be in same directory as `Snakemake` file):

```
snakemake --cores <cores>
```
