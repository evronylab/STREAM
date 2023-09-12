#Concordance statistics for Mendelian trios

The script [trio_concordance.R](scripts/trio_concordance.R) outputs a text file containing genotyping and concordance statistics about the analyzed Mendelian trio called [output.basename]_genotypes_[trioID]_genotypes_concordance_stats.tsv. The statistics calculated in this script are defined here.

The statistics are calculated separately for poly-A and loci with all other motifs (designated "nonA" in these statistics) because poly-A repeats behave differently than loci with longer motifs. Many of the statistics for "nonA" loci are also calculated for each motif length.

Some statistics are calculated separately for HipSTR, GangSTR, and ExpansionHunter to help with assessing the effectiveness of filtering thresholds.


| **STATISTIC** | **DESCRIPTION** | **Also calculated for different motif lengths?** | **Also calculated for each caller?** |
|------------|-----------------|
| num_loci_panel | The number of loci in the panel. |&check;|&cross;|
| num_loci_passing_filters | The number of loci passing all filtering thresholds. |&check;|&cross;|
| num_loci_failing_filters | The number of loci failing filters. |&check;|&cross;|
| num_loci_caller_hipstr | The number of loci whose genotypes are taken from HipSTR. |&cross;|&cross;|
| num_loci_caller_gangstr | The number of loci whose genotypes are taken from GangSTR. |&cross;|&cross;|
| num_loci_caller_eh | The number of loci whose genotypes are taken from ExpansionHunter. |&cross;|&cross;|
| num_loci_fully_genotyped | The number of loci with passing genotypes in all 3 trio members. |&check;|&check;|
| num_concordant_loci | The number of fully-genotyped loci that follow patterns of Mendelian inheritance. |&check;|&check;|
| num_discordant_loci | The number of fully-genotyped loci that do not follow patterns of Mendelian inheritance. |&check;|&check;|
| num_loci_passing_filter_missing_gt | The number of loci that passed filtering but are missing a genotype in one of the trio members. |&check;|&check;|
| frac_num_passing_filters_over_num_loci_panel | The fraction of the panel that passed filtering. |&cross;|&cross;|
| frac_loci_failing_filters_over_num_loci_panel | The fraction of the panel that failed filtering. |&cross;|&cross;|
| frac_num_fully_genotyped_over_num_loci_panel | The fraction of the panel that is fully genotyped. |&cross;|&cross;|
| frac_num_fully_genotyped_over_num_passing_filters | The fraction of loci passing filters that are also fully genotyped. |&cross;|&cross;|
| frac_num_caller_hipstr_over_num_passing_filters | The fraction of loci passing filters whose genotypes are taken from HipSTR. |&cross;|&cross;|
| frac_num_caller_gangstr_over_num_passing_filters | The fraction of loci passing filters whose genotypes are taken from GangSTR. |&cross;|&cross;|
| frac_num_caller_eh_over_num_passing_filters | The fraction of loci passing filters whose genotypes are taken from ExpansionHunter. |&cross;|&cross;|
| frac_loci_passing_filter_missing_gt_over_num_passing_filters | The fraction of loci passing filters that are missing a genotype in one or more trio members. |&check;|&cross;|
| frac_discordant_over_sum_concordant_discordant | The fraction of loci with a valid concordance call that are discordant. |&check;|&check;|