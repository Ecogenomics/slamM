# SDMAss

### Step down metagenome assembler


##### n.b. For ACE users older servers will not work.


### REQUIREMENTS
snakemake

miniconda


### INSTALLATION

#### install miniconda3
```
$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ bash Miniconda3-latest-Linux-x86_64.sh
```


#### install snakemake
```
$ conda install -c bioconda -c conda-forge snakemake
```

#### get the repository

```
$ cd <where_you_want_it_installed>
$ git clone https://github.com/Ecogenomics/SDMass.git
```




### USAGE

The easiest way to use snakemake is to pass configuration options straight from the shell

```
$ cd <install_directory>/SDMass
$ snakemake assemble_pools --cores 24 --use-conda --config long_reads=nanopore_reads.fastq short_reads_1=illumina.1.fastq short_reads_2=illumina.2.fastq workdir=/path/to/output/files/ 
```

If only _short_reads_1_ is given SDMass will assume that the fastq is interleaved paired-end reads

Alternatively you can edit config.yaml to point to reads and directory and then run


```
$ snakemake assemble_pools --cores 12 --use-conda 
```



### KNOWN ISSUES

Unicycler somtimes freezes 

if you see output hang at

```
Aligning reads (2019-04-09 07:20:24)

32,766 / 32,767 (100.0%)
```

n.b. 32,767 will be the number of reads in your particular run. the number preceding (32,766 in this case) will always be one less than this number.

Just exit out and then rerun (Don't worry, it'll pick up where it left off.)