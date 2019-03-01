# SDMAss

### Step down metagenome assembler


### n.b. For ACE users older servers will not work.


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


### USAGE
edit config.yaml to point to reads and directory

cd into SDMass directory

```
$ snakemake step_down_meta_assembly --cores 12 --use-conda 
```

Alternatively files can be passed to snakemake straight from the shell

```
snakemake step_down_meta_assembly --cores 24 --use-conda --config long_reads=reads.fastq workdir=/path/to/output/files/ 
```