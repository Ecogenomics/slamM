# SDMAss

### Step down metagenome assembler


### REQUIREMENTS: 
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
$ snakemake --cores 12 --use-conda step_down_meta_assembly
```