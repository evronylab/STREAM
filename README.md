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
  - [Sequencing data](#sequencing-data)
  - [Samples list](#samples-list)
  - [Formatted lists of microsatellites to analyze](#formatted-lists-of-microsatellites-to-analyze)
  - [List of loci to exclude from ExpansionHunter analysis](#list-of-loci-to-exclude-from-expansionhunter-analysis)
  - [Optional: List of probes used in hybridization capture](#optional-list-of-probes-used-in-hybridization-capture)
  - [Optional: List of supplementary non-microsatellite loci](#optional-list-of-supplementary-non_microsatellite-loci)
- [Pipeline configuration](#pipeline-configuration)
  - [Nextflow configuration file](#nextflow-configuration-file)
  - [Filtering configuration file](#filtering-configuration-file)
- [Running STREAM](#running-pipeline)
  - [Genotype calling](#genotype-calling)
  - [Filtering](#filtering)
  - [Optional: Mendelian discordance](#optional-mendelian-discordance)
- [Outputs](#outputs)
  - [Genotype calling outputs](#genotype-calling-outputs)
  - [Filtering outputs](#filtering-outputs)
  - [Optional: Mendelian discordance outputs](#optional-mendelian-discordance-outputs)
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
 ```[sampleID] [sex] [CRAM] [CRAM index] [trioID] [sampleType]```
 - [sampleID], [sex], [trioID], [sampleType]: Per above
 - [CRAM / CRAM index]: Full path to the CRAM file and its index. The CRAM file must have optical duplicates removed (i.e. removed, not just marked) and be coordinate-sorted.

When given CRAM files as an input, STREAM skips all alignment steps, including optical duplicate removal and sorting. 

### C. Formatted lists of microsatellites to analyze
Several tools in STREAM require a list of microsatellites to analyze. Each tool requires the list in a specific format. We provide [scripts](scripts) for converting a CSV file obtained from [STRATIFY](https://github.com/evronylab/usatShinyApp) into the formats required for the pipeline.

To prepare the microsatellite list files for STREAM, perform the following steps:
1. Prepare a comma-separate file with the fields: row.number,seqnames,start,end,width,period.size,motif,motif.family. When using our [STRATIFY](https://github.com/evronylab/STRATIFY) to design a panel of microsatellites for profiling, this can be downloaded via [STRATIFY's](https://github.com/evronylab/STRATIFY) "View Data" panel. The CSV file will download with a header, but it should be **REMOVED** prior to running the formatting scripts.
   - row.number: A unique ID for each microsatellite locus (string, numeric, or a combination such as MS-#).
   - seqnames: The chromosome on which the microsatellite is located.
   - start/end: The coordinates of the first and last bases in the microsatellite. Coordinates are 1-start, fully closed.
   - width: The length of the microsatellite in base-pairs in the hg38 reference genome.
   - period.size: The length of the microsatellite's repeat unit in base-pairs according to [Tandem Repeats Finder](https://tandem.bu.edu/trf/trf.html).
   - motif: The sequence of the microsatellite's repeat unit.
   - motif.family: The group to which the microsatellite belongs when all possible reverse complements and circular permutations of the motif are reduced to a single motif. E.g. The ATG motif belongs to the ACT motif family because one possible circular permutation of ACT = CAT and the reverse complement of CAT = ATG.

2. Run the following scripts:
   - [convert_to_hipstr.sh](scripts/convert_to_hipstr.sh)
     ```convert_to_hipstr.sh [CSV of targeted microsatellites] [full path to HipSTR-formatted list of microsatellites]```
   - [convert_to_gangstr.sh](scripts/convert_to_gangstr.sh)
     ```convert_to_gangstr.sh [CSV of targeted microsatellites] [full path to GangSTR-formatted list of microsatellites]```
   - [convert_to_eh.sh](scripts/convert_to_eh.sh)
     - Also requires [eh_to_json.R](scripts/eh_to_json.R) to be accessible to Nextflow. eh_to_json.R will be run automatically by convert_to_eh.sh.
     - Requires a list of targeted loci that cause ExpansionHunter to crash and should be [excluded from the ExpansionHunter analysis](#list-of-loci-to-exclude-from-expansionhunter-analysis)
     ```convert_to_eh.sh [CSV of targeted microsatellites] [list of ExpansionHunter exclusion loci] [full path to eh_to_json.R] [full path to ExpansionHunter-formatted list of microsatellites]```
   - [convert_to_bed.sh](scripts/convert_to_bed.sh)
     ```convert_to_bed.sh [CSV of targeted microsatellites] [full path to bedtools-formatted list of microsatellites]```

### D. List of loci to exclude from ExpansionHunter analysis
ExpansionHunter does not process loci with too many Ns in the flanking regions. When this happens, ExpansionHunter stops running instead of skipping the locus. Therefore, we have created a [list](panels/eh_exclusion_list.txt) of loci to exclude from the ExpansionHunter analysis to avoid the analysis ending prematurely. The IDs in the list correspond to the row.number field of the panel CSV when downloaded from [STRATIFY](https://github.com/evronylab/STRATIFY), with the string "MS-" appended as a prefix.

### E. Optional: List of probes used in hybridization capture
If the microsatellites were profiled by hybridization capture, the Picard CollectHsMetrics tool for calculating capture performance metrics requires the coordinates of the targeted microsatellites and the capture probes. The lists for Picard can be prepared as follows:

  - Prepare tab-separated lists of microsatellites and probes with the fields: `[chr] [start coordinate] [end coordinate]`. Coordinates must be 0-based, half-open.
    - A CSV from STRATIFY can be converted with the `awk -F',' -e '{print $2"\t"($3-1)"\t"$4}' list.csv > reformatted_list.txt`
  - Run [convert_to_interval_list.sh](scripts/convert_to_interval_list.sh) on each list
    - The dictionary file format is described in the [Picard BedToIntervalList documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360036883931-BedToIntervalList-Picard)
    ```convert_to_interval_list.sh [microsatellites/probes list] [output file basename] [Reference genome dictionary file]```

Then run [convert_to_interval_list.sh](scripts/convert_to_interval_list.sh)

### F. Optional: List of supplementary non-microsatellite loci
If an additional hybridization capture panel was multiplexed with the microsatellite panel to capture non-microsatellite loci (i.e., a supplementary capture panel), prepare a list of the targeted loci as follows:
  - Prepare tab-separated lists of the targeted loci and probes with the fields: `[chr] [start coordinate] [end coordinate]`. Coordinates must be 0-based, half-open.
  - Run [convert_to_interval_list.sh](scripts/convert_to_interval_list.sh) on each list
    ```convert_to_interval_list.sh [targeted loci/probes list] [output file basename] [Reference genome dictionary file]```

Then run [convert_to_interval_list.sh](scripts/convert_to_interval_list.sh)

## Pipeline configuration
The pipeline is configured with two files:

### A. Nextflow configuration file
A Nextflow configuration file (see [example](config_templates/nextflow.config)) specifying:
- The process executor (i.e. SLURM)
- The input type (params.input_type = 'fastq' or 'cram')
- The 'Samples list' file (params.samples)
- Paths to the reference genome and index
- Paths to input files for tools and scripts
- Output path and file names
- Paths to tools and scripts
- Other runtime parameters

### B. Filtering configuration file
The filtering step reads filtering threshold values from a YAML configuration file (see [example](config_templates/example_filtering_config.yaml)). This configuration file includes:
- Sample-level filters, which are applied to each sample individually
- Locus-level filters, which are evaluated across all samples.
- An order of preference for the three genotype callers. This order determines which caller's genotypes will be used in a situation where the genotype calls from more than one caller pass the filtering thresholds. The second- and third-choice callers can be left blank if desired.
- If performing optional Mendelian concordance analysis, this file can also be used to define the sample names for the father, mother, and child along with the sex of the child.

Filtering thresholds can be set separately for "A" loci (poly-A motif) and "nonA" loci (all other motifs) because poly-A loci are more mutable than loci with longer motifs and may need to be filtered differently.

Note: filters related to read depth can have different values specified for each [sampleType] group as defined in the [samples list](#samples-list), which can be useful when combining different types of samples such as bulk and capture or exome and genome.

## Running STREAM

### A. Genotype calling
The genotype calling step is run with the following command:
```nextflow -C nextflow.config run usat_calling.nf```

### B. Filtering
The genotype calling step outputs an R data frame (RDS file) containing genotype information from all three callers for each sample and locus. The filtering step generates final genotypes by filtering genotype calls on a variety of parameter thresholds and then chooses a caller for each locus whose genotypes are used across all samples at that locus. This is run with the following command:
```Rscript usat_genotype_filtering.R [unfiltered genotypes RDS] [YAML config file]```

### C. Optional: Mendelian discordance
When profiling Mendelian trios (father, mother, child), the filtered RDS file can be analyzed by the [trio_concordance.R](scripts/trio_concordance.R) script to output the Mendelian discordance rate and other relevant statistics. This step is run with the following command:
```Rscript trio_concordance.r [filtered genotypes RDS] [YAML config file]```

## Outputs
Note: [items in brackets] refer to parameters defined in the Nextflow configuration file.

### A. Genotype calling
In the main results directory:
1. nextflow.config: the configuration file used to run STREAM
2. Samples list: the sample list TSV defined in the nextflow.config file
3. [output.basename].rds: R data frame containing genotype information from all three callers for each sample and locus. Columns are defined [here](column_definitions.md).
4. [output.basename]_hipstr.vcf.gz and [output.basename]_hipstr.vcf.gz.csi: output VCF from HipSTR for all samples

In the hipstr_aln_viz directory:
1. [output.basename]_#.aln.viz.gz: graphic rendering of HipSTR's local realignment at each locus. Instructions on how to view the renderings can be found in the [HipSTR documentation](https://hipstr-tool.github.io/HipSTR/#alignment-visualization).

In the sample-specific results directory:
1. [sampleID].cram and [sampleID].cram.crai: aligned CRAM file and index with optical duplicates removed and coordinate-sorted
2. [sampleID].R1_fastqc.html and [sampleID].R2_fastqc.html: fastqc reports (FASTQ start only)
3. [sampleID].markdup.metrics.txt: duplication metrics report, output of Picard MarkDuplicates (FASTQ start only)
4. Picard CollectHsMetrics reports:
   - [sampleID]_picard_usats.txt: output of CollectHsMetrics when the microsatellite coordinates are input as TARGETS and the probe coordinates are input as BAITS
   - [sampleID]_picard_probes.txt: output of CollectHsMetrics when the probe coordinates are input as both TARGETS and BAITS
5. [sampleID]_hipstr.vcf.gz, [sampleID]_gangstr.vcf.gz, [sampleID]_eh.vcf.gz: output VCFs from the three callers
6. [sampleID]_bedtools.txt: output of bedtools coverage; the only columns retained are the name of the locus and the number of reads overlapping the full interval
7. [sampleID].GATK.g.vcf.gz: output of GATK HaplotypeCaller for supplementary panel
8. [sampleID].GATK.vcf.gz: output of GATK GenotypeGVCFs for supplementary panel

### B. Filtering
1. YAML config file: the configuration file defining the filtering thresholds used by the script
2. [output.basename]_genotypes.rds: same data frame as [output.basename].rds with additional columns indicating if the call passed filtering, which filters failed (if any), the caller whose genotypes were used, and the final chosen genotype call. Columns are defined [here](column_definitions.md).

### C. Mendelian discordance
1. [output.basename]\_genotypes_[trioID]_concordance.rds: same data frame as [output.basename]_genotypes.rds with an additional column indicating if the call follows patterns of Mendelian inheritance. Only the samples from the analyzed trio are present in this data frame. Columns are defined [here](concordance_stats_definitions.md).
2. [output.basename]\_genotypes_[trioID]_concordance_stats.tsv: text file listing genotyping and Mendelian concordance statistics for the analyzed trio.

## Citation
If you use STREAM, please cite:

Loh CA\*, Shields DA*, Schwing A, Evrony GD. High-fidelity, large-scale targeted profiling of microsatellites. bioRxiv 2023.month.day.#####; doi: LINK (2023).

Willems, T., Zielinski, D., Yuan, J., Gordon, A., Gymrek, M., & Erlich, Y. (2017). Genome-wide profiling of heritable and de novo STR variations. Nature Methods, 14(6), 590-592. doi:[10.1038/nmeth.4267](https://doi.org/10.1038/nmeth.4267)

Mousavi, N., Shleizer-Burko, S., Yanicky, R., & Gymrek, M. (2019). Profiling the genome-wide landscape of tandem repeat expansions. Nucleic Acids Research, 47(15), e90. doi:[10.1093/nar/gkz501](https://doi.org/10.1093/nar/gkz501)

Dolzhenko, E., Deshpande, V., Schlesinger, F., Krusche, P., Petrovski, R., Chen, S., . . . Eberle, M. A. (2019). ExpansionHunter: a sequence-graph-based tool to analyze variation in short tandem repeat regions. Bioinformatics, 35(22), 4754-4756. doi:[10.1093/bioinformatics/btz431](https://doi.org/10.1093/bioinformatics/btz431)
