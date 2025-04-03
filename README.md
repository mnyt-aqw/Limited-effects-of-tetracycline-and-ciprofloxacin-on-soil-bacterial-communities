# Limited-effects-of-tetracycline-and-ciprofloxacin-on-soil-bacterial-communities

## Read based analysis

The entire read based analysis was performed in a Nextflow pipeline found in the `read_based_pipeline` directory.  

A Nextflow pipeline for analyzing antimicrobial resistance genes (ARGs) in soil metagenomes, with focus on their association with mobile genetic elements.

### Gene Detection

We screen the metagenomes for ARGs using Diamond and a combination of hte CARD and ResFinderFG databases. We also screen for point mutations using a tool called MuMaMe together with the sns.txt file from CARD. Finally we also screen the metagenomes for MGEs using Diamond and a filtered version of the mobileOG-db. You can see how we filtered the database in this notebook `filter_MGE.ipynb`.  

### Taxonomic analysis

In order to determine the taxonomic composition of each sample, as well as extracting the number of 16s sequences we used a tool called metaxaQR.

### Determine if an ARG was mobile

WE downloaded ResSeq plasmid and genomic islands from Island viewer 4 and screened them for ARGs using Diamond. IF a gene was detected it was classified as *evidence for mobility*.

## Genome Resolved Metagenomics

This document describes the bioinformatics workflow used to analyze the effects of tetracycline and ciprofloxacin on soil bacterial communities through genome-resolved metagenomics.

### Assembly and Mapping

First, metagenomic reads from each sample series were co-assembled using Megahit version 1.2.9, which is designed to handle large and complex metagenomes efficiently.

```bash
megahit --presets meta-large \
-1 {SAMPLE}_1_R1_val_1.fq.gz,{SAMPLE}_2_R1_val_1.fq.gz,{SAMPLE}_3_R1_val_1.fq.gz,{SAMPLE}_4_R1_val_1.fq.gz,{SAMPLE}_5_R1_val_1.fq.gz \
-2 {SAMPLE}_1_R2_val_2.fq.gz,{SAMPLE}_2_R2_val_2.fq.gz,{SAMPLE}_3_R2_val_2.fq.gz,{SAMPLE}_4_R2_val_2.fq.gz,{SAMPLE}_5_R2_val_2.fq.gz \
-o {SAMPLE} 
```

The `--presets meta-large` parameter optimizes the assembly process for metagenomic data with high species complexity, adjusting k-mer sizes and other parameters automatically.

After assembly, Bowtie2 was used to map all reads back to their respective assemblies to calculate coverage information, which is essential for binning:

```bash
bowtie2 -x {index} -1 {R1_Reads} -2 {R2_Reads} -S {SAMPLE}.sam --threads 12
samtools view -@ 12 -b {SAMPLE}.sam > {SAMPLE}.bam
samtools sort -@ 12 {SAMPLE}.bam -o {SAMPLE}.sorted.bam
```

The alignment files were converted from SAM to BAM format and sorted for efficient processing in subsequent steps.

### Generate MAGs (Metagenome-Assembled Genomes)

The depth file required for binning was generated using JGI's utility tool:

```bash
jgi_summarize_bam_contig_depths --outputDepth {sample}_depth.txt *.bam
```

This tool calculates the coverage of each contig across all samples, which is crucial information for accurately binning contigs into genomes.

MetaBAT2 was then used to bin contigs into draft genomes (MAGs) based on tetranucleotide frequency and coverage information:

```bash
metabat2 -i {ASSEMBLY} -a {SAMPLE}_depth.txt -o {SAMPLE} --numThreads 12
```

To remove redundancy across samples, dRep was used to dereplicate the MAGs based on genome similarity:

```bash
dRep dereplicate dir_out/ -g *.fa -comp 50 -con 10 -p 10 
```

The parameters `-comp 50 -con 10` set quality thresholds requiring MAGs to have at least 50% completeness and less than 10% contamination to be included in the final set.

### Analyze MAGs

To prepare for multi-sample analysis, a scaffold-to-bin (STB) file was created to track which contigs belong to which MAGs:

```bash
# Generate file with all MAGs
cat *.fa > allGenomes.fasta
parse_stb.py --reverse -f *.fa -o allGenomes.stb
```

Prior to this step, sample names were prepended to contig names to prevent naming collisions between assemblies.

Next, reads from all metagenomes were mapped to the consolidated MAG set:

```bash
bowtie2-build allGenomes.fasta allGenomes
bowtie2 -p 12 -x allGenomes -1 {R1_Reads} -2 {R2_Reads} > {SAMPLE}.sam
samtools view --threads 12 -S -b {SAMPLE}.sam > {SAMPLE}.bam
samtools sort -@ 12 {SAMPLE}.bam -o {SAMPLE}_sorted.bam
```

inStrain's profile module was used to calculate abundance and genetic diversity metrics for each MAG across all samples:

```bash
inStrain profile {SAMPLE}.bam allGenomes.fasta -o {SAMPLE} -p 5 -s allGenomes.stb --database_mode
```

The `--database_mode` parameter optimizes the analysis for large reference collections, which is appropriate for our set of MAGs.

### Taxonomic Classification

GTDB-Tk was used for standardized taxonomic classification of the MAGs according to the Genome Taxonomy Database:

```bash
# Download the latest GTDB database
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_v2_data.tar.gz 
tar xvzf gtdbtk_v2_data.tar.gz

# Run GTDB-Tk classification workflow
gtdbtk classify_wf --genome_dir {INPUT_DIR} --out_dir {OUTPUT_DIR} --cpus 15 --extension .fa
```

GTDB-Tk places MAGs into a taxonomic framework based on concatenated protein phylogeny and ANI comparisons to reference genomes.

### Antibiotic Resistance Gene Detection

To identify potential antibiotic resistance genes (ARGs) in our MAGs:

```bash
# Predict coding regions
prodigal -i {MAG} -a {MAG_NAME}_protein.translations.faa -p single

# Screen for ARGs using CARD database
diamond makedb --in ${CARD} -d CARD_db
diamond blastp --db CARD_db.dmnd --query {MAG_NAME}_protein.translations.faa --out {MAG_NAME}.tsv --threads 1 --id 90 --subject-cover 80 --query-cover 80 --outfmt 6 -k0
```

Prodigal was used to predict protein-coding sequences in each MAG. These predicted proteins were then searched against the Comprehensive Antibiotic Resistance Database (CARD) using DIAMOND, a high-performance BLAST alternative. The search parameters (`--id 90 --subject-cover 80 --query-cover 80`) were set to be relatively stringent to minimize false positives, requiring 90% sequence identity and 80% coverage of both query and subject sequences.

## Data analysis

The entire data analysis can be found in in this notebook `data_analysis.ipynb`. To recreate the environment you can build a Docker container using the `Dockerfile`.