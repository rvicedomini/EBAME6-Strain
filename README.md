# Strain-level metagenome assembly with Strainberry

## 0. Preliminaries

For this tutorial it is recommended to deploy an __Ubuntu 20.04__ VM, possibly with 8 CPUs and at least 16 GB of RAM.

After you connected to your VM through SSH run the following command:
```bash
tmux
```
which will start a new session of Tmux (terminal multiplexer).

Tmux is a handy tool for running one or more terminals at the same time, 
with the possiblity to detach them (so that they keep running in the background) and/or reattach them to a different terminal.

In this tutorial, we will simply use it to prevent assemblers to stop running in case the connection to the VM "breaks".
In the event of this happening, simply reconnect via SSH to your VM and run the command
```bash
tmux attach
```
in order to restore your previous session.


### Download and unpack tutorial data

The two following commands will download the data archive and extract its content in the directory `EBAME6-Strain`:
```bash
wget https://dl.dropbox.com/s/rs56p3sezk84r1o/EBAME6-Strain.tar.gz
tar -xf EBAME6-Strain.tar.gz
```

### Install auxiliary tools

Install some necessary tools in your base conda environment
```bash
conda activate base
conda install -y flye==2.9 purge_dups==1.2.5 minimap2==2.22 samtools==1.14 mummer==3.23 checkm-genome==1.1.3 biopython==1.78
```

### Installing Strainberry

To install Strainberry in the home directory of your VM, run the following commands:
```bash
cd ${HOME}
git clone https://github.com/rvicedomini/strainberry.git
conda env create -n sberry --file ./strainberry/environment.yml
```
which will download Strainberry source code in the `strainberry` directory and set up a dedicated conda environment for Strainberry's dependencies.

To make the `strainberry` command available (without specifying the path to the executable) you can run the following command:
```bash
export PATH=${HOME}/strainberry:${PATH}
```
The previous command should be run each time you open a terminal (a good practice would be to put it in the `~/.bashrc` file so that it would be automatically run each time you open a bash terminal)

## 1. Build Strainberry input

As presented during the lecture, Strainberry requires a strain-oblivious assembly and a long-read alignment as input.
In this section you will assemble a small (mock) metagenome based on _Staphylococcus aureus_.

Before starting, go to the directory containing the data for this tutorial:
```bash
cd ${HOME}/EBAME6-Strain
```

### Use metaFlye to create a strain-oblivious reference assembly (~10 minutes)

First, make sure to activate the conda environment in which you installed the auxiliary tools (the _base_ one in this tutorial).
Then, the simplest way to run metaFlye with a set of long reads in FASTQ/FASTA format is with the following command:
```bash
flye --meta --pacbio-raw ./fastq/saureus_reads.fastq.gz --out-dir ./assemblies/metaflye --threads 8
```
where 
- `./fastq/saureus_reads.fastq.gz` is the path to the reads in FASTQ format
- `./assemblies/metaflye` is the directory where metaFlye will store all output files


In about 10 minutes metaFlye should be able to assemble the reads provided in input. You will focus on two output files:
- `assembly.fasta`: the assembled metagenome in FASTA format
- `assembly_graph.gfa`: the assembly graph in GFA format

<details>
<summary>In the event you did not manage to run metaFlye successfully...</summary>
you can download an assembly of metaFlye for the same dataset along with the corresponding assembly graph:

```bash
# metaFlye assembly
wget https://dl.dropbox.com/s/hfg2vq3biqvphpw/metaflye.fa
# metaFlye assembly graph
wget https://dl.dropbox.com/s/uhf6jwjmjb19ayi/metaflye_graph.gfa
```
</details>


#### Purge duplicated sequences (optional)

metaFlye is a good metagenome assembler that is able to provide high-quality species-level assemblies.
However, in certain cases, it tries to resolve (assemble) strain-specific sequences.
To make the metaFlye assembly as strain-oblivious as possible for Strainberry, 
it is possible to use the tool `purge_dups` to identify/discard some "duplicated" sequences for which metaFlye did not produce a single consensus.

```bash
cd ${HOME}/EBAME6-Strain/assemblies/metaflye
minimap2 -x map-pb assembly.fasta ${HOME}/EBAME6-Strain/fastq/saureus_reads.fastq.gz -t 8 > read_alignment.paf
pbcstat read_alignment.paf
calcuts PB.stat > cutoffs
split_fa assembly.fasta > split.fa
minimap2 -x asm5 -DP split.fa split.fa > self-aln.paf
purge_dups -2 -T cutoffs -c PB.base.cov self-aln.paf > dups.bed
get_seqs dups.bed assembly.fasta
cd ${HOME}/EBAME6-Strain
```
The metaFlye assembly in which duplicated sequences have been purged can be then found in the `${HOME}/EBAME6-Strain/assemblies/metaflye/purged.fa` file.
If you completed this optional step, you should use `purged.fa` as input for the following steps (and not the metaFlye `assembly.fasta` file).


### Align long reads against the strain-oblivious assembly

We are now going to create the second input required by Strainberry, which is a long-read mapping.
We will use the same reads used to generate metaFlye assembly.

```bash
minimap2 -a -x map-pb -t 4 \
  ${HOME}/EBAME6-Strain/assemblies/metaflye/assembly.fasta \
  ${HOME}/EBAME6-Strain/fastq/saureus_reads.fastq.gz \
    | samtools sort --threads 4 -o ${HOME}/EBAME6-Strain/alignments/metaflye_alignment.bam
```
The `minimap2` command will produce an alignment file in the SAM format (`-a` option), with an input consisting of PacBio reads (`-x map-pb` option), and using 4 threads (`-t 4` option).
The `samtools sort` command will sort the alignments by coordinate while outputting a BAM file.


## 2. Run Strainberry to separate strains


Strainberry requires input reference and BAM files to be indexed:
```bash
samtools faidx ${HOME}/EBAME6-Strain/assemblies/metaflye/assembly.fasta
samtools index ${HOME}/EBAME6-Strain/alignments/metaflye_alignment.bam
```

Make sure to activate Strainberry's conda environment:
```bash
conda activate sberry
```

Now you can run Strainberry with input file you generated in the previous step:
```bash
strainberry -r ${HOME}/EBAME6-Strain/assemblies/metaflye/assembly.fasta \
  -b ${HOME}/EBAME6-Strain/alignments/metaflye_alignment.bam \
  -o ${HOME}/EBAME6-Strain/assemblies/strainberry \
  -c CPUs
```
where `CPUs` is the number of CPUs to use (set it according to the virtual machine you deployed for this tutorial)

Once Strainberry finished, it is possible to deactivate the corresponding conda environment:
```bash
conda deactivate
```

<details>
<summary>In the event you did not manage to run Strainberry successfully...</summary>
<br>
you can download a possible output of Strainberry:

```bash
wget https://dl.dropbox.com/s/iqzmgd80iw92jo6/strainberry.fa
```
</details>


## 3. While Strainberry is running... (it could take 15-20 minutes)

On your local laptop, download and install the two following tools:
- [Bandage](https://rrwick.github.io/Bandage/): a program for visualizing _de novo_ assembly graphs
- [Integrative Genomics Viewer (IGV)](https://software.broadinstitute.org/software/igv/download): a program for visualizing an alignment (in BAM format) against a reference sequence and identify, _e.g._, single-nucleotide variants, contig coverage, etc.

Have a look at the assembly graph of metaFlye (file `assembly_graph.gfa` in metaFlye output directory) using Bandage.
The following command will download the assembly graph file in your current (local) directory
(replace `134.158.xxx.xxx` with the actual address of the VM you deployed for this tutorial session):
```bash
scp ubuntu@134.158.xxx.xxx:/home/ubuntu/EBAME6-Strain/assemblies/metaflye/assembly_graph.gfa ./
```
Do you think there are multiple strains of the same species in the dataset?

Now use IGV to look at the long-read alignment against the metaFlye "reference" assembly (i.e., Strainberry's BAM and FASTA input files).
More precisely, load the metaFlye assembly clicking on _Genomes/Load Genome from File..._ and the alignment track clicking on _File/Load from File..._.
Have a look at some of the assembled sequences. Can you identify single-nucleotide variants? Do you think there are multiple strains in this dataset?


## 4. Analysis of the results

### Map strain-oblivious and strain-separated assemblies to a known reference

Use `minimap2` align the assembly of metaFlye and Strainberry to one of the _S. aureus_ genomes in the `references` directory, using the following command:
```bash
minimap2 -ax asm20 REFERENCE.fasta ASSEMBLY.fasta | samtools sort -o OUTPUT.bam
```
where `REFERENCE.fasta`, `ASSEMBLY.fasta`, and `OUTPUT.bam` should be replaced with the paths to the reference, assembly, and output files respectively.

Visualize and compare the two alignments using IGV (it is possible to load multiple BAM files for the same reference).
What do you observe?


### Assembly evaluation metrics

In the `script` directory there is a python script that uses the tool MUMmer to evaluate a metagenome assembly with respect to a given set of reference genomes.
You can run it in the following way:
```bash
python3 script/assembly_stats.py -f ASSEMBLY.fasta -r REFERENCE_DIR -o OUTPUT_DIR
```
where:
- `ASSEMBLY.fasta` is the path to the metagenome assembly to be evaluated (in FASTA format)
- `REFERENCE_DIR` is the path to the directory containing the reference genomes in FASTA format (and with `.fasta` or `.fa` extension)
- `OUTPUT_DIR` is the path to the directory to store the comparison between the assembly and the reference

Run the evaluation script on metaFlye and Strainberry assemblies, using the two _S. aureus_ strain references that are in the directory `${HOME}/EBAME6-Strain/references`.
The script uses MUMmer to partition the sequences of the input assembly according to their best-mapping reference.
More precisely, in the output directory there will be the following files:
- `assembly.saureusATCC.fa`: sequences of `ASSEMBLY.fasta` that have a better match to strain `saureusATCC`
- `assembly.saureusFDAA.fa`: sequences of `ASSEMBLY.fasta` that have a better match to strain `saureusFDAA`
- `assembly.none.fa`: sequences of `ASSEMBLY.fasta` that have no matches
- `report.tsv`: TSV file containing a MUMmer-based evaluation of the assembled sequences

### CheckM analysis (if there is time)

Put metaFlye and Strainberry assemblies in a directory along with the "partitioned" version of Strainberry assembly (obtained from the evaluation in the previous step).
Make sure all files have the `.fa` extension, and run checkm on that directory:
```bash
checkm lineage_wf -x fa -t 8 --tab_table -f assemblies.checkm.tsv --reduced_tree assemblies checkm_output
```


