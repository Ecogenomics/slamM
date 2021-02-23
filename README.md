# slamM

### Step-down metagenome assembler
 
This tool is no-longer being actively developed, though minor updates/bug-fixes may be made.
For ACE users, it is recommended you use newer servers.

### Requirements
- snakemake (last tested with ver. 5.32.0)
- miniconda (last tested with ver. 1.1)


### Installation

#### Install miniconda3 and snakemake (if required)
```
$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ bash Miniconda3-latest-Linux-x86_64.sh
$ conda install -c bioconda -c conda-forge snakemake
```

#### Download the repository

```
$ cd <where_you_want_it_installed>
$ git clone https://github.com/Ecogenomics/slamM.git
```


### Usage

The easiest way to use snakemake is to pass configuration options straight from the shell

```
$ cd <install_directory>/slamm
$ snakemake create_webpage --cores 24 --use-conda --config long_reads=nanopore_reads.fastq short_reads_1=illumina.1.fastq short_reads_2=illumina.2.fastq workdir=/path/to/output/files/
 
```

If only _short_reads_1_ is given slamM will assume that the fastq is interleaved paired-end reads

Alternatively you can edit config.yaml to point to reads and directory and then run


```
$ snakemake create_webpage --cores 24 --use-conda 
```