# microsatellite_calling

PIPELINE is a computational workflow for calling microsatellite genotypes from next-generation sequencing data using three callers: HipSTR, GangSTR, and ExpansionHunter. This repository contains the calling pipeline, quality filtering scripts, and instructions for creating the input files required by each caller.

The calling pipeline can start from either FASTQ or CRAM files. The pipeline also requires a list of microsatellite loci to analyze.


### Outline
- [Computing environment] [#computing-environment]
- [Formatting list of microsatellites] (#formatting-list-of-microsatellites)
- [Running PIPELINE] (#running-pipeline)
  - [Script requirements] (#script-requirements)
  - [Pipeline configuration] (#pipeline-configuration)
  - [Template scripts] (#template-scripts)
- [Quality filtering]
- [Citation] (#citation)

## Computing environment
The pipeline as presented in this manual is set up to be used on a SLURM computing cluster. However, the pipeline can be adapted to work with other job management services that Nextflow supports such as AWS, Azure, and Google Cloud. See the [Nextflow documentation] (https://www.nextflow.io/docs/latest/executor.html#) for more information.


## Formatting list of microsatellites

Several tools in PIPELINE require a list of microsatellites to analyze. Each tool requires the list in a specific format. Here we provide [scripts] (LINK) for converting a CSV file obtained from [APP] (https://github.com/evronylab/usatShinyApp) into the formats required for the pipeline. The input files must be generated prior to running PIPELINE.

The Picard toolkit's CollectHsMetrics also requires the coordinates of the probes used to capture the targeted microsatellites. The list of probes must be in tab-separated format with the following fields: ```[chr] [start coordinate] [end coordinate] [name]```. Coordinates must be 1-start, fully closed.

1. Download CSV list of targeted microsatellites from APP with the following fields: row.number,seqnames,start,end,width,period.size,motif,motif.family
  - The row.number field is used to give unique names to the microsatellites in the list. It can be modified to contain an additional string (i.e., adding "MS-" in front of the row number) or replaced with an your preferred naming system.
2. Run the following scripts:
  - [convert_to_hipstr.sh] (LINK)
  - [convert_to_gangstr.sh] (LINK)
  - [convert_to_eh.sh] (LINK)
  - [convert_to_bed.sh] (LINK)
  - [convert_to_interval_list.sh] (LINK)
3. If there is a supplementary capture panel, run [supplementary_panel_to_interval_list.sh] (LINK)
  - The supplementary panel file must be a TSV file with fields: ```[chr] [start coordinate] [end coordinate]```. Coordinates must be 0-based, half-open.


## Running PIPELINE

### Script requirements
- [fastqc] (https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc)
- [bwa] (https://github.com/lh3/bwa)
- [samtools] (http://www.htslib.org/download/)
- [GATK] (https://github.com/broadinstitute/gatk/releases)
- [Picard .jar file] (https://github.com/broadinstitute/picard/releases/tag/3.0.0)
- [HipSTR] (https://github.com/HipSTR-Tool/HipSTR/releases/tag/v0.6.2)
- [GangSTR] (https://github.com/gymreklab/GangSTR/releases/tag/v2.5)
- [ExpansionHunter] (https://github.com/Illumina/ExpansionHunter/releases/tag/v5.0.0)
- [bedtools] (https://github.com/arq5x/bedtools2/releases/tag/v2.31.0)
- [bcftools] (http://www.htslib.org/download/)
- [R] (https://cran.r-project.org/bin/windows/base/)

### Pipeline configuration
The pipeline can be started from either paired FASTQ or CRAM files. In the config file, params.input_type must be set to 'fastq' or 'cram' to specify your initial data format.

If starting from FASTQ files, your samples must be listed in a TSV file with the following fields (see [example] (LINK)):
  ```[sampleID] [sex] [read 1 FASTQ] [read 2 FASTQ] [trioID] [sampleType] [opticalDistance]```
  - [trioID]: character string specifying the name of the Mendelian trio to which the sample belongs
  - [sampleType]: distinguishes samples from each other by type (i.e., blood vs. skin, exome vs. whole genome) for filtering purposes. Can be "typeA" or "typeB".
  - [opticalDistance]: preferred setting for Picard MarkDuplicates option OPTICAL_DUPLICATE_PIXEL_DISTANCE, which depends on the sequencer used to generate the data.

Starting from CRAM files will cause the pipeline to skip all alignment processes, including removing optical duplicates and sorting. Therefore, any CRAM file used as input for PIPELINE should have optical duplicates removed and be coordinate sorted.

If starting from CRAM files, your samples must be listed in a TSV file with the following fields (see [example] (LINK)):
  ```[sampleID] [sex] [CRAM] [CRAM index] [trioID] [sampleType]```
  
### Template scripts
The processes JOIN_CALLS and JOIN_SAMPLES both call external R scripts. The location of these scripts must be specified in the config file as params.script_path.

## Quality filtering
PIPELINE outputs an R data frame containing genotype information from all three callers for each sample and locus (the data frame is stored as an RDS file). In order to generate a list of final genotypes, the raw data is passed through a script that filters the genotype calls on a variety of quality metrics, then chooses a caller whose genotypes will be used across all samples at each locus.

The filtering script reads filtering threshold values from a YAML configuration file (see [example] (LINK)), then checks if each row of the data frame passed these thresholds. There are sample-level filters, which are applied to each sample individually, and locus-level filters, which are evaluated across all samples. Filters related to read depth can have different values specified for sample type A vs B. 

The user will also need to specify in the YAML an order of preference for the three genotype callers; this order determines which caller's genotypes will be used in a situation where the genotype calls from more than one caller pass the filtering thresholds. The second- and third-choice callers can be left blank if desired.

## Outputs
Note: [items in brackets] refer to parameters defined in the nextflow config file.

### PIPELINE
In the main results directory:
1. [output.basename].rds: R data frame containing genotype information from all three callers for each sample and locus.
2. [output.basename]_genotypes.rds: same data frame as [output.basename].rds with additional columns indicating if the call passed filtering, which filters failed (if any), the caller whose genotypes were used, and the final chosen genotype call.

In the sample-specific results directory:
1. 



## Citation
If you use PIPELINE, please cite:

Loh CA, Shields DA, Schwing A, Truong TK, Evrony GD. High-fidelity, large-scale targeted profiling of microsatellites. bioRxiv 2023.month.day.#####; doi: LINK (2023).

