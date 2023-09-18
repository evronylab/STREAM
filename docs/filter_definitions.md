# STREAM quality filters
The raw genotype data from the genotype calling step must be filtered for quality. Filters are applied to each caller's genotypes separately to account for differences in accuracy and differences in the statistics reported in the VCF files. The tables below defines the filters in the YAML config file used in the filtering script (see [example](https://github.com/evronylab/STREAM/blob/main/config_templates/example_filtering_config.yaml)).

Each of the below filters in the filter YAML config file is defined separately for "A" loci (poly-A motif) with prefix "A." and for "nonA" loci (all other motifs) with prefix "nonA.". This is because poly-A loci are more mutable than loci with longer motifs and may need to be filtered differently.

Additionally, the following filters that are related to read depth can be specified for each [sampleType] group (typeA or typeB) as defined in the Samples list TSV file. This can be useful when performing an analysis combining different types of samples such as non-capture and capture samples, or exome and genome samples:
 - HipSTR: min.mallreads, min.allele.mallreads, min.allreads, min.allele.allreads, min.gldiff, min.mean.mallreads, min.mean.allreads
 - GangSTR: min.total.reads, min.allele.reads, min.mean.total.reads, min.total.reads, min.allele.reads, min.mean.total.reads
 - EH (ExpansionHunter): min.total.reads, min.allele.reads, min.mean.total.reads, min.total.reads, min.allele.reads, min.mean.total.reads

## HipSTR filters:

| **FILTER** | **DEFINITION** | **ASSOCIATED COLUMN(S)** | **FILTER TYPE** |
|------|-------|--------|
| min.qual | Minimum quality score. | H_Q | sample-level |
| min.mallreads | Minimum number of informative reads used to make the genotype call under the maximum likelihood model. | H_MALLREADS_SUM | sample-level |
| min.allele.mallreads | Minimum number of informative reads supporting each allele under the maximum likelihood model. Both alleles must pass the filtering threshold. | H_MALLREADS_1/2 | sample-level |
| min.mvaf | Minimum variant allele fraction (VAF) of each allele under the maximum likelihood model. Both alleles must pass the filtering threshold. | H_MVAF_1/2 | sample-level |
| min.allreads | Minimum number of informative reads used to make the genotype call in the Needleman-Wunsch alignment. | H_ALLREADS_SUM | sample-level |
| min.allele.allreads | Minimum number of informative reads supporting each allele in the Needleman-Wunsch alignment. Both alleles must pass the filtering threshold. | H_ALLREADS_1/2 | sample-level |
| min.avaf | Minimum variant allele fraction (VAF) of each allele in the Needleman-Wunsch alignment. Both alleles must pass the filtering threshold. | H_AVAF_1/2 | sample-level |
| max.stutter.frac | Maximum fraction of reads containing stutter. | H_STUTTERFRAC | sample-level |
| max.flankindel.frac | Maximum fraction of reads containing indels in the flanking sequences. | H_FLINDELFRAC | sample-level |
| min.gldiff | Minimum genotype likelihood difference for each call. | H_GLDIFF | sample-level |
| min.mean.mallreads | Minimum average number of informative reads used to make the genotype call under the maximum likelihood model across all samples. | mean_mdepth_H | locus-level |
| min.mean.allreads | Minimum average number of informative reads used to make the genotype call in the Needleman-Wunsch alignment across all samples. | mean_adepth_H | locus-level |
| min.mean.qual | Minimum average quality score across all samples. | mean_qual_H | locus-level |
| min.sample.frac | Minimum fraction of samples that pass sample-level HipSTR filters. | sample_frac_filter_H | locus-level |

## GangSTR filters:

| **FILTER** | **DEFINITION** | **ASSOCIATED COLUMN(S)** | **FILTER TYPE** |
|------|-------|--------|
| min.qual | Minimum quality score. | G_Q | sample-level |
| min.total.reads | Minimum number of spanning reads at the locus. | G_ENCLREADS_SUM | sample-level |
| min.allele.reads | Minimum number of spanning reads supporting each allele. Both alleles must pass the filtering threshold. | G_ENCLREADS_1/2 | sample-level |
| min.vaf | Minimum variant allele fraction (VAF) of each allele. Both alleles must pass the filtering threshold. | G_VAF_1/2 | sample-level |
| min.mean.total.reads | Minimum average number of spanning reads at the locus across all samples. | mean_depth_G | locus-level |
| min.mean.qual | Minimum average quality score across all samples. | mean_qual_G | locus-level |
| min.sample.frac | Minimum fraction of samples that pass sample-level GangSTR filters. | sample_frac_filter_G | locus-level |

## ExpansionHunter filters:

| **FILTER** | **DEFINITION** | **ASSOCIATED COLUMN(S)** | **FILTER TYPE** |
|------|-------|--------|
| min.total.reads | Minimum number of spanning reads at the locus as measured using bedtools. | B_depth | sample-level |
| min.allele.reads | Minimum number of spanning reads supporting each allele. Both alleles must pass the filtering threshold. | E_ADSP_1/2 | sample-level |
| min.vaf | Minimum variant allele fraction (VAF) of each allele. Both alleles must pass the filtering threshold. | E_VAF_1/2 | sample-level |
| max.vaf | Maximum variant allele fraction (VAF) of each allele. May be > 1 because B_depth is used as the denominator. Both alleles must pass the filtering threshold. | E_VAF_1/2 | sample-level |
| min.lc | Minimum average locus coverage. | E_LC | sample-level |
| min.mean.total.reads | Minimum average number of spanning reads at the locus across all samples. | mean_depth_E | locus-level |
| min.sample.frac | Minimum fraction of samples that pass sample-level ExpansionHunter filters. | sample_frac_filter_E | locus-level |

