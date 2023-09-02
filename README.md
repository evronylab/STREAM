# microsatellite_calling

PIPELINE is a computational workflow for ensemble genotyping of microsatellites from next-generation sequencing data using three callers: HipSTR, GangSTR, and ExpansionHunter. This repository contains instructions for creating the required input files, the genotype calling pipeline, and quality filtering scripts.

The basic steps of PIPELINE from sequencing data to filtered genotype calls are:
1. Format list of microsatellites to be used as input for the various tools in PIPELINE
2. Analyze the sequencing data to obtain unfiltered genotypes from all callers
3. Filter the output and pick final genotypes for each locus

### Outline
- [Computing environment](#computing-environment)
- [Formatting list of microsatellites](#formatting-list-of-microsatellites)
- [Running PIPELINE](#running-pipeline)
  - [Script requirements](#script-requirements)
  - [Pipeline configuration](#pipeline-configuration)
  - [External scripts](#external-scripts)
- [Quality filtering](#quality-filtering)
- [Outputs](#outputs)
- [Citation](#citation)

## Computing environment
The pipeline runs via **Nextflow [link]** and in this manual is set up to be used on a SLURM computing cluster. However, the pipeline can be adapted to work with other job management services that Nextflow supports such as AWS, Azure, and Google Cloud. See the [Nextflow documentation](https://www.nextflow.io/docs/latest/executor.html#) for more information.

## Required tools
The following tools must be available. The paths to them is then specified in the Nextflow configuration file.
- [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc)
- [bwa](https://github.com/lh3/bwa)
- [samtools](http://www.htslib.org/download/)
- [GATK](https://github.com/broadinstitute/gatk/releases)
- [Picard .jar file](https://github.com/broadinstitute/picard/releases/tag/3.0.0)
- [HipSTR](https://github.com/HipSTR-Tool/HipSTR/releases/tag/v0.6.2)
- [GangSTR](https://github.com/gymreklab/GangSTR/releases/tag/v2.5)
- [ExpansionHunter](https://github.com/Illumina/ExpansionHunter/releases/tag/v5.0.0)
- [bedtools](https://github.com/arq5x/bedtools2/releases/tag/v2.31.0)
- [bcftools](http://www.htslib.org/download/)
- [R](https://cran.r-project.org/bin/windows/base/)
- R packages: **

## Input files
### A. Sequencing data
PIPELINE can begin from either FASTQ or CRAM files.

FASTQ files must be **requirements**.

CRAM files must be **requirements**.

It also requires a list of microsatellite loci to analyze, which can be obtained using **APP**--our tool for microsatellite panel design. **When the data is obtained via hybridization capture, PIPLINE can also be given a list of capture probes to calculate hybridization capture performance metrics.**

### B. Samples list
If starting from FASTQ files, prepare a tab-separated file with the following fields (see [example](config_templates/trio_samples_fastq.tsv)):
```[sampleID] [sex] [read 1 FASTQ] [read 2 FASTQ] [trioID] [sampleType] [opticalDistance]```
- [sampleID]: **
- [sex]: **
- [read 1/2 FASTQ]: **
- [trioID]: character string specifying the name of the Mendelian trio to which the sample belongs
- [sampleType]: distinguishes samples from each other by type (i.e., blood vs. skin, exome vs. whole genome) for filtering purposes. Can be "typeA" or "typeB".
- [opticalDistance]: preferred setting for Picard MarkDuplicates option OPTICAL_DUPLICATE_PIXEL_DISTANCE, which depends on the sequencer used to generate the data.

If starting from CRAM files, prepare a tab-separated file with the following fields (see [example](config_templates/trio_samples_cram.tsv)):
 ```[sampleID] [sex] [CRAM] [CRAM index] [trioID] [sampleType]```
 - [sampleID], [sex], [trioID], [sampleType]: per above
 - [CRAM / CRAM index]: full path to the CRAM file and its index. The CRAM file must have optical duplicates removed (i.e. removed, not just marked) and be coordinate-sorted.

When given CRAM files as an input, **PIPELINE** skips all alignment steps, including optical duplicate removal and sorting. 

### C. Formatted lists of microsatellites to analyze
Several tools in PIPELINE require a list of microsatellites to analyze. Each tool requires the list in a specific format. We provide [scripts](scripts) for converting a CSV file obtained from [APP](https://github.com/evronylab/usatShinyApp) into the formats required for the pipeline.

To prepare the microsatellite list files for PIPELINE, perform the following steps:
1. Prepare a comma-separate file with the fields: row.number,seqnames,start,end,width,period.size,motif,motif.family. When using our **APP [link]** to design a panel of microsatellites for profiling, this can be downloaded via **APP**'s **panel**.
   - row.number: A unique ID for each microsatellite locus (string, numeric, or a combination such as MS-#).
   - seqnames:
   - start/end:
   - width:
   - period.size: 
   - motif:
   - motif.family: 

3. Run the following scripts in order:
   - [convert_to_hipstr.sh](scripts/convert_to_hipstr.sh)
   - [convert_to_gangstr.sh](scripts/convert_to_gangstr.sh)
   - [convert_to_eh.sh](scripts/convert_to_eh.sh)
   - [convert_to_bed.sh](scripts/convert_to_bed.sh)
   - [convert_to_interval_list.sh](scripts/convert_to_interval_list.sh)

### D. List of loci to exclude from ExpansionHunter analysis
**Fully explain reason for this, including how the IDs relate to above.

### E. Optional: List of probes used in hybridization capture
If the microsatellites were profiled by hybridization capture, the Picard CollectHsMetrics tool for calculating capture performance metrics requires the coordinates of the capture probes in a tab-separated file with the fields: `[chr] [start coordinate] [end coordinate] [name]`. Coordinates must be 1-start, fully closed. [name] is **.

### F. Optional: List of supplementary non-microsatellite loci
If an additional hybridization capture panel was multiplexed with the microsatellite panel to capture non-microsatellite loci (i.e., a supplementary capture panel), prepare a tab-separated file with the fields: `[chr] [start coordinate] [end coordinate]`. Coordinates must be 0-based, half-open.

Then run [supplementary_panel_to_interval_list.sh](scripts/supplementary_panel_to_interval_list.sh)

## Pipeline configuration
The pipeline is configured with two files:

### A. Nextflow configuration file
A Nextflow configuration file (see [example](nextflow/nextflow.config)) specifying:
- Paths to tools and scripts
- The input type (params.input_type = 'fastq' or 'cram')
- The 'Samples list' file (params.samples)
- Output path and file names
- Other runtime parameters

### B. Filtering configuration file
The filtering step reads filtering threshold values from a YAML configuration file (see [example](config_templates/example_filtering_config.yaml)). This configuration file includes:
- Sample-level filters, which are applied to each sample individually
- Locus-level filters, which are evaluated across all samples.
- An order of preference for the three genotype callers. This order determines which caller's genotypes will be used in a situation where the genotype calls from more than one caller pass the filtering thresholds. The second- and third-choice callers can be left blank if desired.

Note: filters related to read depth can have different values specified for sample type A vs B, which can be useful when combining different types of samples such as bulk and capture or exome and genome.

## Running PIPELINE

### A. Genotype calling
The genotype calling step is run with the following command:
```nextflow -C nextflow.config run usat_calling.nf```

### B. Filtering
The genotype calling step outputs an R data frame (RDS file) containing genotype information from all three callers for each sample and locus. The filtering step generates final genotypes by filtering genotype calls on a variety of parameter thresholds and then chooses a caller for each locus whose genotypes are used across all samples at that locus.

### C. Optional: Mendelian discordance
When profiling Mendelian trios (father, mother, child), the filtered RDS file can be analyzed by the[trio_concordance.R](scripts/trio_concordance.R) script to output the Mendelian discordance rate and other relevant statistics. This step is run with the following command:
**

## Outputs
Note: [items in brackets] refer to parameters defined in the Nextflow configuration file.

### A. Genotype calling
In the main results directory:
1. [output.basename].rds: R data frame containing genotype information from all three callers for each sample and locus.
2. [output.basename]_hipstr.vcf.gz and [output.basename]_hipstr.vcf.gz.csi: output VCF from HipSTR with all samples

In the sample-specific results directory:
1. [sampleID].cram and [sampleID].cram.crai: aligned CRAM file and index with optical duplicates removed and coordinate-sorted
2. [sampleID].R1_fastqc.html and [sampleID].R2_fastqc.html: fastqc reports (FASTQ start only)
3. [sampleID].markdup.metrics.txt: duplication metrics report, output of Picard MarkDuplicates (FASTQ start only)
4. Picard CollectHsMetrics reports:
   - [sampleID]_picard_usats.txt: output of CollectHsMetrics when the microsatellite coordinates are input as TARGETS and the probe coordinates are input as BAITS
   - [sampleID]_picard_probes.txt: output of CollectHsMetrics when the probe coordinates are input as both TARGETS and BAITS
5. [sampleID]_hipstr.vcf.gz, [sampleID]_gangstr.vcf.gz, [sampleID]_eh.vcf.gz: output VCFs from the three callers
6. [sampleID]_bedtools.txt: output of bedtools coverage; the only columns retained are the name of the locus and the number of reads overlapping the full interval
7. [sampleID]_hipstr_table.txt, [sampleID]_gangstr_table.txt, [sampleID]_eh_table.txt: output of BCFtools query for each caller
8. [sampleID].GATK.g.vcf.gz: output of GATK HaplotypeCaller for supplementary panel
9. [sampleID].GATK.vcf.gz: output of GATK GenotypeGVCFs for supplementary panel

### B. Filtering
1. [output.basename]_genotypes.rds: same data frame as [output.basename].rds with additional columns indicating if the call passed filtering, which filters failed (if any), the caller whose genotypes were used, and the final chosen genotype call.
2. [output.basename]_genotypes_[trioID]_genotypes_concordance.rds: same data frame as [output.basename]_genotypes.rds with additional columns indicating if the call follows patterns of Mendelian inheritance. Only the samples from the analyzed trio are present in this data frame.
3. [output.basename]_genotypes_[trioID]_genotypes_concordance_stats.tsv: text file listing genotyping and Mendelian concordance statistics for the analyzed trio.

### C. Mendelian discordance

## Citation
If you use PIPELINE, please cite:

Loh CA, Shields DA, Schwing A, Evrony GD. High-fidelity, large-scale targeted profiling of microsatellites. bioRxiv 2023.month.day.#####; doi: LINK (2023).
