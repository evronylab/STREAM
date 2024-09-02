# STREAM

**STREAM** (**S**hort **T**andem **R**epeat **E**nsemble **A**nalysis **M**ethod) is a computational workflow for ensemble genotyping of microsatellites from next-generation sequencing data using three callers: [HipSTR](https://hipstr-tool.github.io/HipSTR/), [GangSTR](https://github.com/gymreklab/GangSTR), and [ExpansionHunter](https://github.com/Illumina/ExpansionHunter). This repository contains instructions for creating the required input files, the genotype calling pipeline, and quality filtering scripts.

The basic steps of STREAM from sequencing data to filtered genotype calls are:
1. Format list of microsatellites to be used as input for the various tools in STREAM
2. Analyze the sequencing data to obtain unfiltered genotypes from all callers
3. Filter the output and pick final genotypes for each locus

### Outline
- [Computing environment](#computing-environment)
- [Required tools](#required-tools)
- [Input files](#input-files)
  - [Sequencing data](#a-sequencing-data)
  - [Samples list](#b-samples-list)
  - [Formatted lists of microsatellites to analyze](#c-formatted-lists-of-microsatellites-to-analyze)
  - [Optional files if microsatellites were profiled by hybridization capture](#d-optional-files-if-microsatellites-were-profiled-by-hybridization-capture)
  - [Optional files for profiling supplementary non-microsatellite loci](#e-optional-files-for-profiling-supplementary-non-microsatellite-loci)
- [Pipeline configuration](#pipeline-configuration)
  - [Nextflow configuration file](#a-nextflow-configuration-file)
  - [Filtering configuration file](#b-filtering-configuration-file)
- [Running STREAM](#running-stream)
  - [Genotype calling](#a-genotype-calling)
  - [Filtering](#b-filtering)
  - [Optional: Mendelian discordance](#c-optional-mendelian-discordance)
- [Outputs](#outputs)
  - [Genotype calling outputs](#a-genotype-calling-outputs)
  - [Filtering outputs](#b-filtering-outputs)
  - [Mendelian discordance outputs](#c-mendelian-discordance-outputs)
- [Citation](#citation)

## Computing environment
STREAM runs via [Nextflow](https://www.nextflow.io/) and in this manual is set up to be used on a SLURM computing cluster. However, STREAM can be adapted to work with other job management services that Nextflow supports such as AWS, Azure, and Google Cloud. See the [Nextflow documentation](https://www.nextflow.io/docs/latest/executor.html#) for more information.

## Required tools
The following tools must be available. The paths to them is then specified in the Nextflow configuration file.
- [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc)
- [bwa](https://github.com/lh3/bwa)
- [samtools](http://www.htslib.org/download/)
- [GATK](https://github.com/broadinstitute/gatk/releases)
- [Picard .jar file](https://github.com/broadinstitute/picard/releases)
- [HipSTR](https://github.com/HipSTR-Tool/HipSTR/releases)
- [GangSTR](https://github.com/gymreklab/GangSTR/releases)
- [ExpansionHunter](https://github.com/Illumina/ExpansionHunter/releases)
- [bedtools](https://github.com/arq5x/bedtools2/releases)
- [bcftools](http://www.htslib.org/download/)
- [R](https://cran.r-project.org)
- [cutadapt](https://cutadapt.readthedocs.io) (If performing transcriptome filtering, see below)
- [STAR](https://github.com/alexdobin/STAR) (If performing transcriptome filtering, see below)
- [seqkit](https://bioinf.shenwei.me/seqkit) (If performing transcriptome filtering, see below)
- R packages:
  - [dplyr](https://cloud.r-project.org/web/packages/dplyr/index.html)
  - [tidyr](https://cloud.r-project.org/web/packages/tidyr/index.html)
  - [stringr](https://cloud.r-project.org/web/packages/stringr/index.html)
  - [configr](https://cran.r-project.org/web/packages/configr/index.html)

## Input files
### A. Sequencing data
STREAM can begin from either FASTQ or CRAM files.

FASTQ files must be paired and have read names that follow the regular expression pattern `[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9\-]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*`.

CRAM files must be coordinate-sorted with optical duplicates removed. Duplicate marking should also be removed from the CRAM before using in the pipeline.

### B. Samples list
If starting from FASTQ files, prepare a tab-separated file with the following fields (see [example](config_templates/trio_samples_fastq.tsv)):
```[sampleID] [sex] [read 1 FASTQ] [read 2 FASTQ] [trioID] [sampleType] [opticalDistance]```
- [sampleID]: The unique ID given to the sample.
- [sex]: The sex of the sample. Must be listed as either 'M' or 'F'.
- [read 1/2 FASTQ]: Full paths to the paired FASTQ files for the sample. If there is > 1 FASTQ per sample, the FASTQs can be written in a comma-separated list; the read 1 and 2 FASTQs must be written in the same order to maintain proper pairing. 
- [trioID]: Character string specifying the name of the Mendelian trio to which the sample belongs
- [sampleType]: Distinguishes samples from each other by type (i.e., blood vs. skin, exome vs. whole genome) for filtering purposes. Can be "typeA" or "typeB".
- [opticalDistance]: Setting for Picard MarkDuplicates option OPTICAL_DUPLICATE_PIXEL_DISTANCE, which depends on the sequencer used to generate the data.

If starting from CRAM files, prepare a tab-separated file with the following fields (see [example](config_templates/trio_samples_cram.tsv)):
 `[sampleID] [sex] [CRAM] [CRAM index] [trioID] [sampleType]`
- [sampleID], [sex], [trioID], [sampleType]: Per above
- [CRAM / CRAM index]: Full path to the CRAM file and its index. The CRAM file must have optical duplicates removed (i.e. removed, not just marked) and be coordinate-sorted.

When given CRAM files as an input, STREAM skips all alignment steps, including optical duplicate removal and sorting. 

### C. Formatted lists of microsatellites to analyze
STREAM requires files listing the microsatellite loci to analyze. Each tool in STREAM requires the list of loci in a specific format. We provide [scripts](scripts) for converting a main microsatellite metadata CSV file obtained from [STRATIFY](https://github.com/evronylab/STRATIFY) into the various file formats required for the pipeline.

To prepare all the required formatted lists of microsatellites for STREAM, perform the following steps:

1. Prepare the main microsatellite metadata CSV file with the fields: `row.number,seqnames,start,end,width,period.size,motif,motif.family`.
   - This file is specified for STREAM in the params.panel field of the Nextflow configuration file
   - When using our [STRATIFY](https://github.com/evronylab/STRATIFY) tool to design a panel of microsatellites for profiling, this can be downloaded via STRATIFY's "View Data" panel. The CSV file will download with a header, but it should be **REMOVED** prior to running the below formatting scripts.
   - The fields in the CSV file are:
     - row.number: A unique ID for each microsatellite locus (string, numeric, or a combination such as MS-#).
     - seqnames: The chromosome on which the microsatellite is located.
     - start/end: The coordinates of the first and last bases in the microsatellite. Coordinates are 1-start, fully closed.
     - width: The length of the microsatellite in base-pairs in the hg38 reference genome.
     - period.size: The length of the microsatellite's repeat unit in base-pairs according to [Tandem Repeats Finder](https://tandem.bu.edu/trf/trf.html).
     - motif: The sequence of the microsatellite's repeat unit.
     - motif.family: The group to which the microsatellite belongs when all possible reverse complements and circular permutations of the motif are reduced to a single motif. E.g. The ATG motif belongs to the ACT motif family because one possible circular permutation of ACT = CAT and the reverse complement of CAT = ATG.

2. Convert the main microsatellite metadata CSV file to other file formats required by STREAM using the following scripts:

   a. [convert_to_hipstr.sh](scripts/convert_to_hipstr.sh): `convert_to_hipstr.sh [main microsatellites metadata CSV] [output file name of HipSTR-formatted list of microsatellites]`
      - The output file of this script is then specified for STREAM in the params.hipstr_usats field of the Nextflow configuration file

   b. [convert_to_gangstr.sh](scripts/convert_to_gangstr.sh): `convert_to_gangstr.sh [main microsatellites metadata CSV] [output file name of GangSTR-formatted list of microsatellites]`
      - The output file of this script is then specified for STREAM in the params.gangstr_usats field of the Nextflow configuration file

   c. [convert_to_eh.sh](scripts/convert_to_eh.sh): `convert_to_eh.sh [main microsatellites metadata CSV] [list of ExpansionHunter exclusion loci] [full path to eh_to_json.R] [output file name of ExpansionHunter-formatted list of microsatellites]`
      - The output file of this script is then specified for STREAM in the params.eh_usats field of the Nextflow configuration file
      - This script requires the [eh_to_json.R](scripts/eh_to_json.R) script.
      - This script requires a list of loci that cause ExpansionHunter to crash and that should be excluded from the ExpansionHunter analysis. This is because loci with too many Ns in the flanking regions cause ExpansionHunter to stop running instead of skipping the locus. We have created a [list](panels/eh_exclusion_list.txt) of loci to exclude from the ExpansionHunter analysis to avoid the analysis ending prematurely. The IDs in the list correspond to the row.number field of the panel CSV when downloaded from [STRATIFY](https://github.com/evronylab/STRATIFY), with the string "MS-" appended as a prefix.
	
   d. [convert_to_bed.sh](scripts/convert_to_bed.sh): `convert_to_bed.sh [main microsatellites metadata CSV] [output file name of bedtools-formatted list of microsatellites]`
      - The output file of this script is then specified for STREAM in the params.bedtools_usats field of the Nextflow configuration file

### D. Optional files if microsatellites were profiled by hybridization capture
If the microsatellites were profiled by hybridization capture, the Picard CollectHsMetrics tool for calculating capture performance metrics requires two files, one containing the coordinates of the targeted microsatellites and one containing the coordinates of the capture probes. These files are prepared as follows:

a. Make two tab-separated files of coordinates, one of the microsatellite loci and one of the probes, with the fields: `[chr] [start coordinate] [end coordinate]`. Coordinates must be 0-based, half-open (BED format).
   - Note: the file of microsatellite loci in this format can be created from the main microsatellite metadata CSV file obtained from STRATIFY with the command: `awk -F',' -e '{print $2"\t"($3-1)"\t"$4}' [STRATIFY_loci.csv] > [STRATIFY_loci.bed]`

b. Run [convert_to_interval_list.sh](scripts/convert_to_interval_list.sh) on each file with the command: `convert_to_interval_list.sh [microsatellites/probes list] [output file basename] [Reference genome dictionary file]`
   - The output is [output file basename].interval_list
   - Note: the dictionary (.dict) file format is described in the [Picard BedToIntervalList documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360036883931-BedToIntervalList-Picard)

The resulting microsatellite and probe coordinate files are then specified for STREAM in the params.usats_panel_picard and params.usats_probes_picard fields, respectively, of the Nextflow configuration file.

### E. Optional files for profiling supplementary non-microsatellite loci
Supplementary non-microsatellite loci can be genotyped with GATK by providing to STREAM coordinates of the supplementary loci. If these supplementary loci were profiled by an additional hybridization capture panel that was multiplexed with the microsatellite panel, capture metrics can be calculated for the supplementary loci panel with the Picard CollectHsMetrics tool by providing STREAM the coordinates of the capture probes.

These files are prepared as follows:

a. Make tab-separated files, one of the supplementary loci and if needed, one of the probes, with the fields: `[chr] [start coordinate] [end coordinate]`. Coordinates must be 0-based, half-open (BED format).

b. Run [convert_to_interval_list.sh](scripts/convert_to_interval_list.sh) on each file with the command: `convert_to_interval_list.sh [targeted loci/probes list] [output file basename] [Reference genome dictionary file]`
   - The output is [output file basename].interval_list
   - Note: the dictionary (.dict) file format is described in the [Picard BedToIntervalList documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360036883931-BedToIntervalList-Picard)

The resulting supplementary loci and probe coordinate files are then specified for STREAM in the params.supp_panel_picard and params.supp_probes_picard fields, respectively, of the Nextflow configuration file. Note, if the supplementary loci were not targeted by hybridization capture, the user can specify only the params.supp_panel_picard file.

## Pipeline configuration
The pipeline is configured with two files:

### A. Nextflow configuration file
A Nextflow configuration file (see [example](config_templates/nextflow.config)) specifying:
- The process executor (i.e. SLURM)
- The input type (params.input_type = 'fastq' or 'cram')
- The 'Samples list' file (params.samples)
- Paths to the reference genome and index
- Paths to reference genome for STAR transcriptome alignment for 'transcriptome filtering' (see below). This can be generated with the command: STAR --runMode genomeGenerate --genomeDir [genomedir] --genomeFastaFiles [genome fasta] --sjdbGTFfile [GENCODE gtf file] --sjdbOverhang 100"
- Paths to loci to analyze
- Output path and file names
- Paths to tools and scripts
- Other runtime parameters

### B. Filtering configuration file
The filtering step reads filtering threshold values from a YAML configuration file (see [example](config_templates/example_filtering_config.yaml)). The parameters defined by the configuration file are documented [here](docs/filter_definitions.md), and include:
- Sample-level filters, which are applied to each sample individually.
- Locus-level filters, which are evaluated across all samples.
- An order of preference for the three genotype callers. This order determines which caller's genotypes will be used in a situation where the genotype calls from more than one caller pass the filtering thresholds. The second- and third-choice callers can be left blank if desired.
- If performing optional Mendelian concordance analysis, the sample names for the father, mother, and child, and the sex of the child.

Filtering thresholds can be set separately for "A" loci (poly-A motif) and "nonA" loci (all other motifs) because poly-A loci are more mutable than loci with longer motifs and may need to be filtered differently.

Note: filters related to read depth can have different values specified for each [sampleType] group (typeA or typeB) as defined in the [samples list](#b-samples-list), which can be useful when performing an analysis combining different types of samples such as non-capture and capture samples, or exome and genome samples.

## Running STREAM

### A. Genotype calling
The genotype calling step applies all three callers to produce initial unfiltered genotype calls for each sample and locus. The main output is an RDS data frame containing the genotype calls.

The genotype calling can optionally perform 'transcriptome filtering' of FASTQ reads prior to alignment for STR genotype calling. This is performed by aligning reads to the transcriptome and then filtering out those reads that align as spliced reads. This can be useful when analyzing DNA sequencing reads from multi-omic RNA+DNA single-cell samples that have some residual RNA-derived reads in the DNA data.

The genotype calling step is run with the following command:
`nextflow -C nextflow.config run usat_calling.nf`

### B. Filtering
The filtering step generates final genotypes by filtering calls from the genotype calling step on a variety of parameter thresholds and then choosing a caller for each locus whose genotypes are used across all samples at that locus. 

The filtering step is run with the following command:
`Rscript usat_genotype_filtering.R [RDS output of the genotype calling step] [YAML config file]`

### C. Optional: Mendelian discordance
When profiling Mendelian trios (father, mother, child), the filtered RDS file can be analyzed by the [trio_concordance.R](scripts/trio_concordance.R) script to output the Mendelian discordance rate and other relevant statistics. This step is run with the following command:
`Rscript trio_concordance.r [filtered genotypes RDS] [YAML config file]`

## Outputs
Note: [items in brackets] refer to parameters defined in the Nextflow configuration file.

### A. Genotype calling outputs
In the main results directory:
1. nextflow.config: the configuration file used to run STREAM, for documentation
2. Samples list: the sample list TSV defined in the nextflow.config file, for documentation
3. [output.basename].rds: R data frame containing genotype information from all three callers for each sample and locus. Columns are defined [here](docs/column_definitions.md).
4. [output.basename]_hipstr.vcf.gz and [output.basename]_hipstr.vcf.gz.csi: output VCF from HipSTR for all samples

In the hipstr_aln_viz directory:
1. [output.basename]_#.aln.viz.gz: graphic rendering of HipSTR's local realignment at each locus. Instructions on how to view the renderings can be found in the [HipSTR documentation](https://hipstr-tool.github.io/HipSTR/#alignment-visualization).

In the sample-specific results directory:
1. [sampleID].cram and [sampleID].cram.crai: aligned CRAM file and index with optical duplicates removed and coordinate-sorted
2. [sampleID].R1_fastqc.html and [sampleID].R2_fastqc.html: fastqc reports (FASTQ start only)
3. [sampleID].STAR.Log.final.out: Output of STAR transcriptome alignment if performing optional transcriptome filtering
4. [sampleID].STAR.splicemetrics.txt: Fraction of total reads that are splice reads per transcriptome alignment, which are filtered from FASTQ files
5. [sampleID].markdup.metrics.txt: duplication metrics report, output of Picard MarkDuplicates (FASTQ start only)
6. Picard CollectHsMetrics reports:
   - [sampleID]_picard_usats.txt: output of CollectHsMetrics when the microsatellite coordinates are input as TARGETS and the probe coordinates are input as BAITS
   - [sampleID]_picard_probes.txt: output of CollectHsMetrics when the probe coordinates are input as both TARGETS and BAITS
7. [sampleID]_hipstr.vcf.gz, [sampleID]_gangstr.vcf.gz, [sampleID]_eh.vcf.gz: output VCFs from the three callers
8. [sampleID]_bedtools.txt: output of bedtools coverage; the only columns retained are the name of the locus and the number of reads overlapping the full interval
9. [sampleID].GATK.g.vcf.gz: output of GATK HaplotypeCaller, if a supplementary panel was specified
10. [sampleID].GATK.vcf.gz: output of GATK GenotypeGVCFs, if a supplementary panel was specified

### B. Filtering outputs
1. YAML config file: a copy of the configuration file defining the filtering thresholds used by the script, for documentation
2. [output.basename]_genotypes.rds: same data frame as [output.basename].rds with additional columns indicating if the call passed filtering, which filters failed (if any), the caller whose genotypes were used, and the final chosen genotype call. Columns are defined [here](docs/column_definitions.md).

### C. Mendelian discordance outputs
If the optional Mendelian discordance step was performed, the following outputs are produced:
1. [output.basename]\_genotypes_[trioID]_concordance.rds: same data frame as [output.basename]_genotypes.rds with an additional column indicating if the call follows patterns of Mendelian inheritance. Only the samples from the analyzed trio are present in this data frame. Columns are defined [here](docs/concordance_stats_definitions.md).
2. [output.basename]\_genotypes_[trioID]_concordance_stats.tsv: text file listing genotyping and Mendelian concordance statistics for the analyzed trio.

## Citation
If you use STREAM, please cite:

STREAM: Loh CA\*, Shields DA*, Schwing A, Evrony GD. High-fidelity, large-scale targeted profiling of microsatellites. bioRxiv; doi: [https://doi.org/10.1101/2023.11.28.569106](https://www.biorxiv.org/content/10.1101/2023.11.28.569106v1) (2023).

HipSTR: Willems, T., Zielinski, D., Yuan, J., Gordon, A., Gymrek, M., & Erlich, Y. (2017). Genome-wide profiling of heritable and de novo STR variations. Nature Methods, 14(6), 590-592. doi:[10.1038/nmeth.4267](https://doi.org/10.1038/nmeth.4267)

GangSTR: Mousavi, N., Shleizer-Burko, S., Yanicky, R., & Gymrek, M. (2019). Profiling the genome-wide landscape of tandem repeat expansions. Nucleic Acids Research, 47(15), e90. doi:[10.1093/nar/gkz501](https://doi.org/10.1093/nar/gkz501)

ExpansionHunter: Dolzhenko, E., Deshpande, V., Schlesinger, F., Krusche, P., Petrovski, R., Chen, S., . . . Eberle, M. A. (2019). ExpansionHunter: a sequence-graph-based tool to analyze variation in short tandem repeat regions. Bioinformatics, 35(22), 4754-4756. doi:[10.1093/bioinformatics/btz431](https://doi.org/10.1093/bioinformatics/btz431)
