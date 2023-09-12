# Data frame columns

The final outputs of STREAM are R data frames containing basic information about each locus along with results from HipSTR, GangSTR, ExpansionHunter, and bedtools. The columns are defined below.

## Locus and sample information
These columns contain basic information about the locus and sample. The locus information is retrieved from the list of microsatellites input into STREAM, and the sample information is extracted from the STREAM sample list. These columns are present in all output data frames.

| **COLUMN** | **DESCRIPTION** |
|------------|-----------------|
| name | The unique name of the locus (string, numeric, or a combination such as MS-#). |
| chr | The chromosome on which the locus is located. |
| start | The coordinate of the first base of the microsatellite in the reference genome. |
| end | The coordinate of the last base of the microsatellite in the reference genome. |
| width | The length in base-pairs of the microsatellite in the reference genome. |
| period.size | The length in base-pairs of the microsatellite repeat unit according to Tandem Repeats Finder. |
| motif | The sequence of the microsatellite's repeat unit. |
| motif.family | The group to which the microsatellite belongs when all possible reverse complements and circular permutations of the motif are reduced to a single motif. E.g. The ATG motif belongs to the ACT motif family because one possible circular permutation of ACT = CAT and the reverse complement of CAT = ATG. |
| sampleID | The unique ID given to the sample. |
| sex | The sex of the sample, listed as either 'M' or 'F'. |
| trio | Character string specifying the name of the Mendelian trio to which the sample belongs. |
| sampleType | Specifies whether the sample belongs to the groups 'typeA' or 'typeB'. This column allows for two sample types to be filtered with different settings. |
| motifLength | The length of the character string in the motif.family column. This is generally the same as period.size but on rare occasion period.size does not match the length of motif.family. |

## Caller columns
The following columns are extracted directly from the VCF outputs of HipSTR, GangSTR, and ExpansionHunter. The name of the column after the "H/G/E_" prefix directly corresponds to the name of the extracted field from the VCF. The definitions of these columns can be found in the documentation for their respective callers. These columns are present in all output data frames.

### HipSTR columns
[(documentation)](https://hipstr-tool.github.io/HipSTR/#file-formats).
- H_GT
- H_GB (also split into separate columns for each allele, H_GB_1/2)
- H_Q
- H_DP
- H_DSTUTTER
- H_DFLANKINDEL
- H_ALLREADS (also split into separate columns for each allele, H_ALLREADS_1/2)
- H_MALLREADS (also split into separate columns for each allele, H_ALLREADS_1/2)
- H_GLDIFF

### ExpansionHunter columns
[(documentation)](https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/RepeatGenotyping_fDG_dtSW.htm).
- E_REF
- E_GT
- E_REPCN (also split into separate columns for each allele, E_REPCN_1/2)
- E_SO
- E_ADSP (also split into separate columns for each allele, E_ADSP_1/2)
- E_ADFL
- E_ADIR
- E_LC

### GangSTR columns
[(documentation)](https://github.com/gymreklab/GangSTR#formats).
- G_REF
- G_GT
- G_Q
- G_REPCN (also split into separate columns for each allele, G_REPCN_1/2)
- G_REPCI
- G_DP
- G_RC
- G_ENCLREADS (also split into separate columns for each allele, G_ENCLREADS_1/2)
- G_FLNKREADS

## PIPELINE-generated columns
The following columns are generated as part of the data processing in **PIPELINE**. These columns are present in all output data frames.

| **COLUMN** | **DESCRIPTION** |
|------------|-----------------|
| B_depth | The number of reads in the CRAM that overlap the microsatellite coordinates. Corresponds to the first output column of bedtools coverage. |
| E_REPDIFF_1/2 | The difference between the number of repeat units in the reference genome and the number of repeat units in each allele, as called by ExpansionHunter. REPDIFF = REPCN - REF |
| G_REPDIFF_1 | The difference between the number of repeat units in the reference genome and the number of repeat units in each allele, as called by GangSTR. REPDIFF = REPCN - REF |
| H_STUTTERFRAC | The fraction of informative reads that contain stutter, as called by HipSTR. STUTTERFRAC = DSTUTTER / DP |
| H_FLINDELFRAC | The fraction of informative reads that contain indels in the flanks, as called by HipSTR. FLINDELFRAC = DFLANKINDEL / DP |
| H_REPDIFF_1/2 | The difference between the number of repeat units in the reference genome and the number of repeat units in each allele, as called by HipSTR. REPDIFF = GB / motifLength |
| H_REPCN_1/2 | The number of repeat units in each allele, as called by HipSTR. REPCN = (width + GB) / motifLength |
| H_MALLREADS_SUM | The total number of reads listed in the HipSTR MALLREADS field. |
| H_ALLREADS_SUM | The total number of reads listed in the HipSTR ALLREADS field. |
| H_MVAF_1/2 | The variant allele fraction of each allele, as called by HipSTR. MVAF = MALLREADS_1/2 / MALLREADS_SUM |
| H_AVAF_1/2 | The variant allele fraction of each allele, as called by HipSTR. AVAF = ALLREADS_1/2 / ALLREADS_SUM |
| E_VAF_1 | The variant allele fraction of each allele, as called by ExpansionHunter. ExpansionHunter does not report a total read depth, so B_depth is used as the denominator, which may lead to VAF values > 1. VAF = ADSP_1/2 / B_depth |
| G_ENCLREADS_SUM | The total number of reads listed in the GangSTR ENCLREADS field. |
| G_VAF_1/2 | The variant allele fraction of each allele, as called by GangSTR. VAF = ENCLREADS_1/2 / ENCLREADS_SUM |


### Filtering-generated columns
The following columns are generated by the filtering script. These columns are present in [output.basename]_genotypes.rds and [output.basename]_genotypes_[trioID]_genotypes_concordance.rds.

| **COLUMN** | **DESCRIPTION** |
|------------|-----------------|
| HipSTR_passsample | Indicates whether the sample passed sample-level filters for the HipSTR fields. |
| GangSTR_passsample | Indicates whether the sample passed sample-level filters for the GangSTR fields. |
| EH_passsample | Indicates whether the sample passed sample-level filters for the ExpansionHunter fields. |
| sample_frac_filter_H | The fraction of samples that passed sample-level filters for the HipSTR fields. |
| sample_frac_filter_G | The fraction of samples that passed sample-level filters for the GangSTR fields. |
| sample_frac_filter_E | The fraction of samples that passed sample-level filters for the ExpansionHunter fields. |
| mean_mdepth_H | The mean H_MALLREADS_SUM across all samples. This field is calculated separately for different sample types. |
| mean_adepth_H | The mean H_ALLREADS_SUM across all samples. This field is calculated separately for different sample types. |
| mean_depth_G | The mean G_ENCLREADS_SUM across all samples. This field is calculated separately for different sample types. |
| mean_depth_E | The mean B_depth across all samples. This field is calculated separately for different sample types. B_depth is used because ExpansionHunter does not report a total read depth. |
| mean_qual_H | The mean H_Q across all samples. This field is calculated separately for different sample types. |
| mean_qual_G | The mean G_Q across all samples. This field is calculated separately for different sample types. |
| locus_means_H | Indicates whether the locus passed the mean depth and quality filters for HipSTR fields. |
| locus_means_G | Indicates whether the locus passed the mean depth and quality filters for GangSTR fields. |
| locus_means_E | Indicates whether the locus passed the mean depth filter for ExpansionHunter fields. |
| HipSTR_passlocus | Indicates whether the sample passed locus-level filters for the HipSTR fields. |
| GangSTR_passlocus | Indicates whether the sample passed locus-level filters for the GangSTR fields. |
| EH_passlocus | Indicates whether the sample passed locus-level filters for the ExpansionHunter fields. |
| caller | The caller whose genotypes will be used across the entire locus. Can be HipSTR, GangSTR, EH, or "no caller". "no caller" occurs when the locus fails to pass filters for all three callers. |
| allele_1/2 | The difference in length, in base-pairs, between the called alleles and the reference genome. The reported alleles are from the caller listed in the "caller" column. allele_2 is NA for loci on chrX or chrY in male samples. |
| allele_1/2_repcn | The total number of repeat units in the called alleles. The reported alleles are from the caller listed in the "caller" column. allele_2 is NA for loci on chrX or chrY in male samples. |
| sample_failed_filters | List of the sample-level filters that were not passed at each locus. |
| locus_failed_filters | List of the locus-level filters that were not passed across all samples. |

### Concordance calculation columns
The following columns are generated when calculating Mendelian concordance. These columns are present in [output.basename]_genotypes_[trioID]_genotypes_concordance.rds.

| **COLUMN** | **DESCRIPTION** |
|------------|-----------------|
| concordance | Indicates whether the called alleles follow Mendelian inheritance patterns. TRUE = concordant, FALSE = discordant, NA = concordance could not be determined (i.e. if one sample was missing a genotype) |

