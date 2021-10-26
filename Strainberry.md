# Strain-level metagenome assembly with Strainberry

## Preliminaries

It is reccommended to create a VM, possibly with 4-8 CPUs and at least 16 GB of RAM.

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

Follow the instructions at [Strainberry repository](https://github.com/rvicedomini/strainberry#installation) in the **Installation** section.
There is **no** need to run it on the example dataset to verify the correct installation.


## Build Strainberry input

As presented during the lecture, Strainberry requires a strain-oblivious assembly and a long-read alignment as input.

### Use metaFlye to create a strain-oblivious reference assembly

#### Purge duplicated sequences [optional]

metaFlye is a good metagenome assembler that is able to provide high-quality species-level assemblies.
However, in certain cases, it tries to resolve (assemble) strain-specific sequences.
To make the metaFlye assembly as strain-oblivious as possible, 
it is possible to use the tool `purge_dups` to identify/discard some "duplicated" sequences for which metaFlye did not produce a single consensus.

```bash

```


### Align long reads against the strain-oblivious assembly


## Run Strainberry to separate strains

Activate Strainberry's conda environment:
```bash
conda activate sberry
```
Run Strainberry with the assembly and the long-read alignments as input:
```bash
```

## While Strainberry is running... (it could take 15-20 minutes)

On your local laptop, download and install the two following tools:
- [Bandage](https://rrwick.github.io/Bandage/): a program for visualising _de novo_ assembly graphs
- [Integrative Genomics Viewer (IGV)](https://software.broadinstitute.org/software/igv/download): a program for visualizing an alignment (in BAM format) against a reference sequence

Then have a look at the assembly graph of metaFlye (file `assembly_graph.gfa` in metaFlye output directory) using Bandage.
How many strains do you think are present in the assembled dataset?
Now use IGV to look at the long-read alignment against the metaFlye "reference" assembly (i.e., Strainberry's BAM and FASTA input files).
What do you observe?


## Evaluation of the results

### Map strain-oblivious and strain-separated assemblies to a known reference

Use `minimap2` align the assembly of metaFlye and Strainberry to one of the _S. aureus_ genomes in the `references` directory, using the following command:
```bash
minimap2 -ax asm20 [reference.fasta] [assembly.fasta] | samtools sort -o [output.bam] -
```
where `[reference.fasta]`, `[assembly.fasta]`, and `[output.bam]` should be replaced with the reference, assembly, and output files respectively.

Visualize and compare the two alignments using IGV.

### Assembly evaluation metrics

In the `script` directory there is a python script that uses the tool MUMmer to evaluate a metagenome assembly with respect to a given set of reference genomes.
You can run it in the following way:
```bash
python3 script/assembly_stats.py -f [assembly.fasta] -r [reference-dir] -o [output-dir]
```
where
- `[assembly.fasta]` is the metagenome assembly to be evaluated
- `[reference-dir]` is the directory containing the reference genomes in FASTA format (and with `.fasta` or `.fa` extension)
- `[output-dir]` is the output directory to store the results

