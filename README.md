[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
<img align="right" src="docs/images/tetris.png" height="100">
# tetris

**TETRIS** is a Nextflow pipeline for DNA mapping and variant calling. It provides a flexible and scalable workflow for processing sequencing data, from raw reads to variant calls.

It trims reads with ([`fastp`](https://github.com/OpenGene/fastp)), aligns with ([`BWA-MEM`](https://bio-bwa.sourceforge.net/)), marks duplicates (optional) with ([`GATK MarkDuplicates`](https://gatk.broadinstitute.org/hc/en-us/articles/21905036102043-MarkDuplicates-Picard)), and calls variants with ([`BCFTOOLS`](https://www.htslib.org/)). Additionally QC stats are computed with ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)), ([`Samtools`](https://www.htslib.org/)) and ([`mosdepth`](https://github.com/brentp/mosdepth)) which is aggregated into a report by ([`MultiQC`](http://multiqc.info/))


<img align="centre" src="docs/images/tetris-metro-map.png" width="700">

The pipeline is split into two subworkflows, mapping and variant calling. The individual subworkflows can be run depending on your use case, or you can run both as a complete pipeline.

## Usage
To run the pipeline with default parameters:

``` bash
# To run complete pipeline
nextflow run main.nf --input samplesheet.csv --reference reference.fasta

# To only running mapping steps
nextflow run main.nf --input samplesheet.csv --reference reference.fasta --mapping_only

# To only running variant calling steps, start from bams
nextflow run main.nf --input samplesheet.csv --reference reference.fasta --calling_only
```


### Input
The pipeline requires two main inputs:
1. A samplesheet CSV files with the follow columns:
- name: Sample name
- seqid: sequencing ID
- seq_type: 'paired' or 'single' end read data type
- fastq_1: path to forward reads
- fastq_2: path to reverse reads (optional for single end data)

or when starting at the variant calling subworkflow:
- name: Sample name
- bam: path to bam file
- bai: path to bam index file

2. A reference genome in FASTA format

### Parameters
| **Parameter**       | **Description**                                       | **Default**      |
|---------------------|-------------------------------------------------------|------------------|
| --input             | Input samplesheet CSV file                            | (required)       |
| --reference         | Reference genome FASTA file                           | (required)       |
| --outdir            | Output directory                                      | results          |
| --split_fastq       | Number of reads to split FASTQ files by               | 0 (no splitting) |
| --skip_fastp        | Skip fastp processing                                 | false            |
| --skip_markdup      | Skip duplicate read marking                           | false            |
| --markdup_tool      | Tool for marking duplicates                           | gatk             |
| --grouped_call      | Run mpileup and call on all BAMs together             | false            |
| --mapping_only      | Perform only read mapping (no variant calling)        | false            |
| --regions           | bed file of genomic regions to call variants in       | null             |
| --constrain_alleles | Constrain alleles (requires indexed target positions) | false            |
| --targets_index     | Indexed set of target positions and alleles           | null             |

## Output
The pipeline outputs, bams, vcfs, multqc and pipeline reports, will be saved in the specified output directory (default: results).


## TODO:

- add in example data (current example data is too large for GitHub)
- ensure compatibility with single end read data
- add more optional samplesheet info to include in read group header
- add second duplicate read mapping tool option
