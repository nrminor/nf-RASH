# nf-RASH: Nextflow Regional ASsembly Helper

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/) [![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/) [![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

### Overview

RASH is a containerized Nextflow pipeline that helps do two, conceptually simple tasks:

1. Extracting genome regions of interest (e.g., the MHC region) from both PacBio HiFi and Oxford Nanopore sequencing reads.
2. Running them through high-accuracy de novo assembly using [Hifiasm](https://github.com/chhylp123/hifiasm).

These tasks are made considerably more difficult by the fact that the FASTQ files used to store these kinds of reads can be *enormous*, often approaching half a terabyte each. By wrapping these tasks in a Nextflow pipeline, it becomes trivial to parallelize tasks and keep track of intermediate files. Nextflow also makes it possible to containerize the software used at each step (RASH supports Docker and Apptainer/Singularity at this time), allowing users to deploy the pipeline on systems ranging from a Desktop PC to a high-performance computing (HPC) cluster. We typically run RASH on a slurm-based HPC cluster and have included our config as an example at [config/hpc.config](config/hpc.config).

### Usage

To get started, all you need installed is [Nextflow](https://nextflow.io/) and either Docker or Apptainer. To run it on your own data, you use a command like this:

```zsh
nextflow run . \
-profile singularity \
-c config/hpc.config \
--pb_fastq inputs/pbhifi.fastq.gz \
--ont_fastq inputs/ont.fastq.gz \
--ref_fasta inputs/ref.fasta.gz \
--desired_regions resources/desired_regions.tsv \
--split_max 500000 \
--cpus 20
```

This command demonstrates a few ways you can configure the pipeline, namely:

- (`nextflow run .` runs the pipeline in the current working directory)
- `-profile singularity` tells RASH to user the Singularity container engine. This can also be set to `apptainer` or `docker`
- `-c config/hpc.config` tells it to pull additional configuration from our example HPC config. We use additional config files to specify platform-specific configurations, e.g., to interface with the Slurm workload manager.
- `--pb_fastq` and `--ont_fastq` specify paths to input files. RASH will check that these files actually exist. Note also that the folder `inputs/` is included in the repo `.gitignore` file.
- `--ref_fasta` (you guessed it) specifies the path to the reference FASTA file.
- `--desired_regions` expects a TSV file that specifies the genome coordinates for where you'd like to assemble. See our example at [`resources/desired_regions.tsv`](resources/desired_regions.tsv). Ultimately, this gets parsed into a `samtools view` expression.
- `--split_max` specifies the maximum number of reads to include in each split FASTQ. This comes up in the first step in the pipeline, which is to split up the enormous input FASTQs into a bunch of smaller chunks, each of which can be mapped to the reference in parallel.
- `--cpus` specifies the number of cores *per task*. So, if your input FASTQ gets split up into 80 FASTQs, the mapping task for each of those 80 files will *each get 20 cores*. This is of course a setting for HPC clusters. If you're running RASH locally, you can leave this argument out, as the default value here is 6.

### Final Note

De novo assembly is a complex process that requires manual review. Don't be rash--review your assemblies!
