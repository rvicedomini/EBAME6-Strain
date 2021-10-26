# Strain deconvolution with Strainberry

## Preliminaries

It is reccommended to create a VM with 4-8 CPUs and at least 16 GB of RAM.

### Download and unpack tutorial data

The two following commands will download the data archive and extract its content in the directory `EBAME6-Strain`:
```bash
wget https://dl.dropbox.com/s/rs56p3sezk84r1o/EBAME6-Strain.tar.gz
tar -xf EBAME6-Strain.tar.gz
```

### Install auxiliary tools

Install some necessary tools in your conda environment
```bash
conda config --add channels bioconda
conda config --add channels conda-forge
conda install flye purge_dups minimap2 samtools mummer
```

### Installing Strainberry

Follow the instructions at [Strainberry repository](https://github.com/rvicedomini/strainberry)


## Build Strainberry input



### Use metaFlye to create a strain-oblivious reference assembly

#### [optional] Purge duplicated sequences

### Compare first assembly and purged one

### Align long reads on the chosen input reference



## Run Strainberry to separate strains

Activate Strainberry's conda environment:
```bash
conda activate sberry
```
Run Strainberry with the assembly and the long-read alignments as input:
```bash
```


## Evaluate results

### Map assembly to a known reference

### Assembly evaluation metrics

### CheckM evaluation
