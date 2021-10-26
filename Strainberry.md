# Strain-level metagenome assembly with Strainberry

## 0. Preliminaries

It is reccommended to create a VM, possibly with 4--8 CPUs and at least 16 GB of RAM.

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

To install Strainberry in the home directory of your VM, run the following commands:
```bash
git clone https://github.com/rvicedomini/strainberry.git
conda env create -n sberry --file environment.yml
```
which will download Strainberry source code in the `strainberry` directory and set up a dedicated conda environment for Strainberry's dependencies.

To make the `strainberry` command available (without specifying the path to the executable) you can run the following command:
```bash
export PATH=~/strainberry:${PATH}
```

## 1. Build Strainberry input

As presented during the lecture, Strainberry requires a strain-oblivious assembly and a long-read alignment as input.
In this section you will assemble a small (mock) metagenome based on _Staphylococcus aureus_.

Before starting, go to the directory containing the data for this tutorial:
```bash
cd ~/EBAME6-Strain
```

### Use metaFlye to create a strain-oblivious reference assembly

The simplest way to run metaFlye with a set of long reads in FASTQ/FASTA format is with the following command:
```bash
flye --meta --pacbio-raw ./fastqs/saureus_reads.fastq.gz --out-dir ./assemblies/metaflye --threads [CPUs]
```
where `[CPUs]` should be replaced with the number of CPUs you want to use (set it according to your deployed VM).
In about 10 minutes metaFlye should be able to assemble the reads provided in input. You will focus on two output files:
- `assembly.fasta`: the assembled metagenome in FASTA format
- `assembly_graph.gfa`: the assembly graph in GFA format


#### Purge duplicated sequences (optional)

metaFlye is a good metagenome assembler that is able to provide high-quality species-level assemblies.
However, in certain cases, it tries to resolve (assemble) strain-specific sequences.
To make the metaFlye assembly as strain-oblivious as possible for Strainberry, 
it is possible to use the tool `purge_dups` to identify/discard some "duplicated" sequences for which metaFlye did not produce a single consensus.

```bash
cd assemblies/metaflye
minimap2 -x map-pb assembly.fasta ../../fastq/saureus_reads.fastq.gz -t 4 > read_alignment.paf
pbcstat read_alignment.paf
calcuts PB.stat > cutoffs
split_fa assembly.fasta > split.fa
minimap2 -x asm5 -DP split.fa split.fa > self-aln.paf
purge_dups -2 -T cutoffs -c PB.base.cov self-aln.paf > dups.bed
get_seqs dups.bed assembly.fasta
```
The metaFlye assembly in which duplicated sequences have been purged can be then found in the `assemblies/metaflye/purged.fa` file.
If you completed this optional step, you should use `purged.fa` as input for the following steps (and not the metaFlye `assembly.fasta` file).


### Align long reads against the strain-oblivious assembly


## 2. Run Strainberry to separate strains

First, activate Strainberry's conda environment:
```bash
conda activate sberry
```
Move to the `assemblies` directory and run Strainberry with input file you generated in the previous step:
```bash
cd ~/EBAME6-Strain/assemblies
strainberry -r ./metaflye/assembly.fasta -b ../alignments/metaflye_alignment.bam -o sberry_metaflye -c [CPUs]
```
where `[CPUs]` is the number of CPUs to use (set it according to the virtual machine you deployed for this tutorial)

Now it is possible to deactivate the conda environment:
```bash
conda deactivate
```

## While Strainberry is running... (it could take 15-20 minutes)

On your local laptop, download and install the two following tools:
- [Bandage](https://rrwick.github.io/Bandage/): a program for visualizing _de novo_ assembly graphs
- [Integrative Genomics Viewer (IGV)](https://software.broadinstitute.org/software/igv/download): a program for visualizing an alignment (in BAM format) against a reference sequence and identify, e.g., single-nucleotide variants, contig coverage, etc.

Have a look at the assembly graph of metaFlye (file `assembly_graph.gfa` in metaFlye output directory) using Bandage.
Do you think there are multiple strains of the same species in the dataset?

Now use IGV to look at the long-read alignment against the metaFlye "reference" assembly (i.e., Strainberry's BAM and FASTA input files).
More precisely, load the metaFlye assembly clicking on _Genomes/Load Genome from File..._ and the alignment track clicking on _File/Load from File..._.
Have a look at some of the assembled sequences. Can you identify single-nucleotide variants? Do you think there are multiple strains in this dataset?


## Evaluation of the results

### Map strain-oblivious and strain-separated assemblies to a known reference

Use `minimap2` align the assembly of metaFlye and Strainberry to one of the _S. aureus_ genomes in the `references` directory, using the following command:
```bash
minimap2 -ax asm20 [reference.fasta] [assembly.fasta] | samtools sort -o [output.bam] -
```
where `[reference.fasta]`, `[assembly.fasta]`, and `[output.bam]` should be replaced with the reference, assembly, and output files respectively.

Visualize and compare the two alignments using IGV (it is possible to load multiple BAM files for the same reference).

### Assembly evaluation metrics

In the `script` directory there is a python script that uses the tool MUMmer to evaluate a metagenome assembly with respect to a given set of reference genomes.
You can run it in the following way:
```bash
python3 script/assembly_stats.py -f [assembly.fasta] -r [reference-dir] -o [output-dir]
```
where you should replace:
- `[assembly.fasta]` with the metagenome assembly to be evaluated
- `[reference-dir]` with the directory containing the reference genomes in FASTA format (and with `.fasta` or `.fa` extension)
- `[output-dir]` with the name of the output directory to store the evaluation results



