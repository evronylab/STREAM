#!/share/apps/r/4.3.1/gcc/bin/Rscript

# usage: Rscript trio_concordance.r [filtered genotypes RDS] [YAML config file]

# Filtered genotypes RDS is output of usat_genotype_filtering.r with genotype data for trios
# Config file is a YAML file with parameter settings for each caller and ensemble calling

# CONCORDANCE LOGIC:
# If the locus has genotypes for all samples and is concordant, concordance == TRUE
# If the locus has genotypes for all samples and is discordant, concordance == FALSE
# If the locus is missing a genotype for any sample, concordance == "missing genotype"
## Father can be missing genotype for chrX concordance in trios with male child without disrupting concordance calculation
## Mother can be missing genotype for chrY concordance in trios with male child without disrupting concordance calculation
# When concordance column is joined back to full list of loci, any locus that did not pass filtering will have concordance == NA

library(dplyr)
library(tidyr)
library(configr)

args <- commandArgs(trailingOnly = T)

# Load genotypes list from RDS
results <- readRDS(args[1])

# Load YAML config file
config <- read.config(file = args[2])

# filter genotypes for trio of interest
results <- results %>% filter(trio == config$sample$trio)


# save sample name as variable
father_1 <- paste0("allele_1_",config$sample$father)
father_2 <- paste0("allele_2_",config$sample$father)
mother_1 <- paste0("allele_1_",config$sample$mother)
mother_2 <- paste0("allele_2_",config$sample$mother)
child_1 <- paste0("allele_1_",config$sample$child)
child_2 <- paste0("allele_2_",config$sample$child)


# Prepare results data frame for concordance calculation
trio <- results %>%
  # filter for loci with calls
  filter(caller %in% c("HipSTR","GangSTR","EH")) %>%
  # keep call data and some metadata (for joining)
  select(c(name, chr, sampleID, allele_1, allele_2)) %>%
  pivot_wider(names_from = sampleID, values_from = c(allele_1, allele_2))


# Determine concordance at each locus
# Concordance logic is dependent on the sex of the child in the trio
if (config$sample$sex.of.child == "M") {
  
  # calculate by row
  trio <- trio %>% rowwise() %>%
    mutate(concordance =
             case_when(
               # if the chromosome is X and call is concordant
               chr == "chrX" &
                 # if no genotypes are missing
                 !anyNA(c(get(mother_1),get(mother_2),get(child_1))) &
                 # is either maternal allele present in the child?
                 (get(mother_1) == get(child_1) | get(mother_2) == get(child_1)) &
                 # is the first child allele present in the mother?
                 # second allele on X chromosome is always NA in males
                 get(child_1) %in% c(get(mother_1),get(mother_2)) ~ TRUE,
               
               # if the chromosome is Y and call is concordant
               chr == "chrY" &
                 # if no genotypes are missing
                 !anyNA(c(get(father_1),get(child_1))) &
                 # is the first paternal allele present in the child?
                 get(father_1) == get(child_1) ~ TRUE,
               
               # if the chromosome is an autosome and call is concordant
               chr != "chrX" & chr != "chrY" &
                 # if no genotypes are missing
                 !anyNA(c(get(father_1),get(father_2),get(mother_1),get(mother_2),get(child_1),get(child_2))) &
                 # is either paternal allele present in the child?
                 (get(father_1) %in% c(get(child_1),get(child_2)) | get(father_2) %in% c(get(child_1),get(child_2))) &
                 # is either maternal allele present in the child?
                 (get(mother_1) %in% c(get(child_1),get(child_2)) | get(mother_2) %in% c(get(child_1),get(child_2))) &
                 # is the first child allele present in at least one parent?
                 get(child_1) %in% c(get(father_1),get(father_2),get(mother_1),get(mother_2)) &
                 # is the second child allele present in at least one parent?
                 get(child_2) %in% c(get(father_1),get(father_2),get(mother_1),get(mother_2)) ~ TRUE,
               
               
               # if the chromosome is X and any maternal allele or child allele is missing
               chr == "chrX" & anyNA(c(get(mother_1),get(mother_2),get(child_1))) ~ NA,
               
               # if the chromosome is Y and the paternal or child allele is missing
               chr == "chrY" & anyNA(c(get(father_1),get(child_1))) ~ NA,
               
               # if the chromosome is an autosome and any allele is missing
               chr != "chrX" & chr != "chrY" & anyNA(c(get(father_1),get(father_2),get(mother_1),get(mother_2),get(child_1),get(child_2))) ~ NA,
               
               # if the call is discordant
               .default = FALSE))
  
} else {
  
  # calculate by row
  trio <- trio %>% rowwise() %>%
    mutate(concordance =
             case_when(
               # if the chromosome is X and the call is concordant
               chr == "chrX" &
                 # if no genotypes are missing
                 !anyNA(c(get(father_1),get(mother_1),get(mother_2),get(child_1),get(child_2))) &
                 # is either paternal allele present in the child?
                 get(father_1) %in% c(get(child_1),get(child_2)) &
                 # is either maternal allele present in the child?
                 (get(mother_1) %in% c(get(child_1),get(child_2)) | get(mother_2) %in% c(get(child_1),get(child_2))) &
                 # is the first child allele present in at least one parent?
                 get(child_1) %in% c(get(father_1),get(mother_1),get(mother_2)) &
                 # is the second child allele present in at least one parent?
                 get(child_2) %in% c(get(father_1),get(mother_1),get(mother_2)) ~ TRUE,
               
               # if the chromosome is an autosome and the call is concordant
               chr != "chrX" & chr != "chrY" &
                 # if no genotypes are missing
                 !anyNA(c(get(father_1),get(father_2),get(mother_1),get(mother_2),get(child_1),get(child_2))) &
                 # is either paternal allele present in the child?
                 (get(father_1) %in% c(get(child_1),get(child_2)) | get(father_2) %in% c(get(child_1),get(child_2))) &
                 # is either maternal allele present in the child?
                 (get(mother_1) %in% c(get(child_1),get(child_2)) | get(mother_2) %in% c(get(child_1),get(child_2))) &
                 # is the first child allele present in at least one parent?
                 get(child_1) %in% c(get(father_1),get(father_2),get(mother_1),get(mother_2)) &
                 # is the second child allele present in at least one parent?
                 get(child_2) %in% c(get(father_1),get(father_2),get(mother_1),get(mother_2)) ~ TRUE,
               
               # if the chromosome is Y, disregard call because chrY concordance isn't possible with a female child
               chr == "chrY" ~ NA,
               
               # if the chromosome is X and either maternal allele, either child allele, or the paternal allele are missing
               chr == "chrX" & anyNA(c(get(father_1),get(mother_1),get(mother_2),get(child_1),get(child_2))) ~ NA,
               
               # if the chromosome is an autosome and any allele is missing
               chr != "chrX" & chr != "chrY" & anyNA(c(get(father_1),get(father_2),get(mother_1),get(mother_2),get(child_1),get(child_2))) ~ NA,
               
               # if the call is discordant,
               .default = FALSE))
}


results <- left_join(results, select(trio,c(name,concordance)), by = c("name"))

nonA_results_stats <- results %>% group_by(name, chr, motifLength, motif.family, caller, concordance) %>% summarise() %>% filter(motif.family != "A")
A_results_stats <- results %>% group_by(name, chr, motifLength, motif.family, caller, concordance) %>% summarise() %>% filter(motif.family == "A")

# Calculate statistics for concordance and filtering data

# non poly-A loci stats
nonA_stats <-
  data.frame(
    nonA_num_loci_panel = nrow(nonA_results_stats),
    nonA_num_2bp_loci_panel = length(which(nonA_results_stats$motifLength == 2)),
    nonA_num_3bp_loci_panel = length(which(nonA_results_stats$motifLength == 3)),
    nonA_num_4bp_loci_panel = length(which(nonA_results_stats$motifLength == 4)),
    nonA_num_loci_passing_filters = length(which(nonA_results_stats$caller != "no call")),
    nonA_num_2bp_loci_passing_filters = length(which(nonA_results_stats$motifLength == 2 & nonA_results_stats$caller != "no call")),
    nonA_num_3bp_loci_passing_filters = length(which(nonA_results_stats$motifLength == 3 & nonA_results_stats$caller != "no call")),
    nonA_num_4bp_loci_passing_filters = length(which(nonA_results_stats$motifLength == 4 & nonA_results_stats$caller != "no call")),
    nonA_num_loci_failing_filters = length(which(nonA_results_stats$caller == "no call")),
    nonA_num_2bp_loci_failing_filters = length(which(nonA_results_stats$motifLength == 2 & nonA_results_stats$caller == "no call")),
    nonA_num_3bp_loci_failing_filters = length(which(nonA_results_stats$motifLength == 3 & nonA_results_stats$caller == "no call")),
    nonA_num_4bp_loci_failing_filters = length(which(nonA_results_stats$motifLength == 4 & nonA_results_stats$caller == "no call")),
    nonA_num_loci_caller_hipstr = length(which(nonA_results_stats$caller == "HipSTR")),
    nonA_num_loci_caller_gangstr = length(which(nonA_results_stats$caller == "GangSTR")),
    nonA_num_loci_caller_eh = length(which(nonA_results_stats$caller == "EH")),
    nonA_num_loci_fully_genotyped = length(which(nonA_results_stats$concordance == TRUE | nonA_results_stats$concordance == FALSE)),
    nonA_num_2bp_loci_fully_genotyped = length(which((nonA_results_stats$concordance == TRUE | nonA_results_stats$concordance == FALSE)
                                                     & nonA_results_stats$motifLength == 2)),
    nonA_num_3bp_loci_fully_genotyped = length(which((nonA_results_stats$concordance == TRUE | nonA_results_stats$concordance == FALSE)
                                                     & nonA_results_stats$motifLength == 3)),
    nonA_num_4bp_loci_fully_genotyped = length(which((nonA_results_stats$concordance == TRUE | nonA_results_stats$concordance == FALSE)
                                                     & nonA_results_stats$motifLength == 4)),
    nonA_num_hipstr_loci_fully_genotyped = length(which((nonA_results_stats$concordance == TRUE | nonA_results_stats$concordance == FALSE)
                                                        & nonA_results_stats$caller == "HipSTR")),
    nonA_num_gangstr_loci_fully_genotyped = length(which((nonA_results_stats$concordance == TRUE | nonA_results_stats$concordance == FALSE)
                                                         & nonA_results_stats$caller == "GangSTR")),
    nonA_num_eh_loci_fully_genotyped = length(which((nonA_results_stats$concordance == TRUE | nonA_results_stats$concordance == FALSE)
                                                    & nonA_results_stats$caller == "EH")),
    nonA_num_concordant_loci = length(which(nonA_results_stats$concordance == TRUE)),
    nonA_num_2bp_concordant_loci = length(which(nonA_results_stats$motifLength == 2 & nonA_results_stats$concordance == TRUE)),
    nonA_num_3bp_concordant_loci = length(which(nonA_results_stats$motifLength == 3 & nonA_results_stats$concordance == TRUE)),
    nonA_num_4bp_concordant_loci = length(which(nonA_results_stats$motifLength == 4 & nonA_results_stats$concordance == TRUE)),
    nonA_num_hipstr_concordant_loci = length(which(nonA_results_stats$caller == "HipSTR" & nonA_results_stats$concordance == TRUE)),
    nonA_num_gangstr_concordant_loci = length(which(nonA_results_stats$caller == "GangSTR" & nonA_results_stats$concordance == TRUE)),
    nonA_num_eh_concordant_loci = length(which(nonA_results_stats$caller == "EH" & nonA_results_stats$concordance == TRUE)),
    nonA_num_discordant_loci = length(which(nonA_results_stats$concordance == FALSE)),
    nonA_num_2bp_discordant_loci = length(which(nonA_results_stats$motifLength == 2 & nonA_results_stats$concordance == FALSE)),
    nonA_num_3bp_discordant_loci = length(which(nonA_results_stats$motifLength == 3 & nonA_results_stats$concordance == FALSE)),
    nonA_num_4bp_discordant_loci = length(which(nonA_results_stats$motifLength == 4 & nonA_results_stats$concordance == FALSE)),
    nonA_num_hipstr_discordant_loci = length(which(nonA_results_stats$caller == "HipSTR" & nonA_results_stats$concordance == FALSE)),
    nonA_num_gangstr_discordant_loci = length(which(nonA_results_stats$caller == "GangSTR" & nonA_results_stats$concordance == FALSE)),
    nonA_num_eh_discordant_loci = length(which(nonA_results_stats$caller == "EH" & nonA_results_stats$concordance == FALSE)),
    nonA_num_loci_passing_filter_missing_gt = length(which(is.na(nonA_results_stats$concordance) & nonA_results_stats$caller != "no call")),
    nonA_num_2bp_loci_passing_filter_missing_gt = length(which(nonA_results_stats$motifLength == 2 & is.na(nonA_results_stats$concordance) & nonA_results_stats$caller != "no call")),
    nonA_num_3bp_loci_passing_filter_missing_gt = length(which(nonA_results_stats$motifLength == 3 & is.na(nonA_results_stats$concordance) & nonA_results_stats$caller != "no call")),
    nonA_num_4bp_loci_passing_filter_missing_gt = length(which(nonA_results_stats$motifLength == 4 & is.na(nonA_results_stats$concordance & nonA_results_stats$caller != "no call"))),
    nonA_num_hipstr_loci_passing_filter_missing_gt = length(which(nonA_results_stats$caller == "HipSTR" & is.na(nonA_results_stats$concordance))),
    nonA_num_gangstr_loci_passing_filter_missing_gt = length(which(nonA_results_stats$caller == "GangSTR" & is.na(nonA_results_stats$concordance))),
    nonA_num_eh_loci_passing_filter_missing_gt = length(which(nonA_results_stats$caller == "EH" & is.na(nonA_results_stats$concordance)))) %>%
  
  mutate(
    nonA_frac_num_passing_filters_over_nonA_num_loci_panel = nonA_num_loci_passing_filters / nonA_num_loci_panel,
    nonA_frac_loci_failing_filters_over_nonA_num_loci_panel = nonA_num_loci_failing_filters / nonA_num_loci_panel,
    nonA_frac_num_fully_genotyped_over_nonA_num_loci_panel = nonA_num_loci_fully_genotyped / nonA_num_loci_panel,
    nonA_frac_num_fully_genotyped_over_nonA_num_passing_filters = nonA_num_loci_fully_genotyped / nonA_num_loci_passing_filters,
    nonA_frac_num_caller_hipstr_over_nonA_num_passing_filters = nonA_num_loci_caller_hipstr / nonA_num_loci_passing_filters,
    nonA_frac_num_caller_gangstr_over_nonA_num_passing_filters = nonA_num_loci_caller_gangstr / nonA_num_loci_passing_filters,
    nonA_frac_num_caller_eh_over_nonA_num_passing_filters = nonA_num_loci_caller_eh / nonA_num_loci_passing_filters,
    nonA_frac_loci_passing_filter_missing_gt_over_nonA_num_passing_filters = nonA_num_loci_passing_filter_missing_gt / nonA_num_loci_passing_filters,
    nonA_frac_2bp_loci_passing_filter_missing_gt_over_nonA_num_2bp_passing_filters = nonA_num_2bp_loci_passing_filter_missing_gt / nonA_num_2bp_loci_passing_filters,
    nonA_frac_3bp_loci_passing_filter_missing_gt_over_nonA_num_3bp_passing_filters = nonA_num_3bp_loci_passing_filter_missing_gt / nonA_num_3bp_loci_passing_filters,
    nonA_frac_4bp_loci_passing_filter_missing_gt_over_nonA_num_4bp_passing_filters = nonA_num_4bp_loci_passing_filter_missing_gt / nonA_num_4bp_loci_passing_filters,
    nonA_frac_discordant_over_nonA_sum_concordant_discordant = nonA_num_discordant_loci / (nonA_num_concordant_loci + nonA_num_discordant_loci),
    nonA_frac_2bp_discordant_over_nonA_sum_2bp_concordant_discordant = 
      nonA_num_2bp_discordant_loci / (nonA_num_2bp_concordant_loci + nonA_num_2bp_discordant_loci),
    nonA_frac_3bp_discordant_over_nonA_sum_3bp_concordant_discordant =
      nonA_num_3bp_discordant_loci / (nonA_num_3bp_concordant_loci + nonA_num_3bp_discordant_loci),
    nonA_frac_4bp_discordant_over_nonA_sum_4bp_concordant_discordant =
      nonA_num_4bp_discordant_loci / (nonA_num_4bp_concordant_loci + nonA_num_4bp_discordant_loci),
    nonA_frac_hipstr_discordant_over_nonA_sum_hipstr_concordant_discordant =
      nonA_num_hipstr_discordant_loci / (nonA_num_hipstr_concordant_loci + nonA_num_hipstr_discordant_loci),
    nonA_frac_gangstr_discordant_over_nonA_sum_gangstr_concordant_discordant =
      nonA_num_gangstr_discordant_loci / (nonA_num_gangstr_concordant_loci + nonA_num_gangstr_discordant_loci),
    nonA_frac_eh_discordant_over_nonA_sum_eh_concordant_discordant =
      nonA_num_eh_discordant_loci / (nonA_num_eh_concordant_loci + nonA_num_eh_discordant_loci))


# poly A loci stats  
A_stats <-
  data.frame(
    A_num_loci_panel = nrow(A_results_stats),
    A_num_loci_passing_filters = length(which(A_results_stats$caller != "no call")),
    A_num_loci_failing_filters = length(which(A_results_stats$caller == "no call")),
    A_num_loci_caller_hipstr = length(which(A_results_stats$caller == "HipSTR" )),
    A_num_loci_caller_gangstr = length(which(A_results_stats$caller == "GangSTR" )),
    A_num_loci_caller_eh = length(which(A_results_stats$caller == "EH" )),
    A_num_loci_fully_genotyped = length(which(A_results_stats$concordance == TRUE | A_results_stats$concordance == FALSE)),
    A_num_hipstr_loci_fully_genotyped = length(which((A_results_stats$concordance == TRUE | A_results_stats$concordance == FALSE)
                                                     & A_results_stats$caller == "HipSTR" )),
    A_num_gangstr_loci_fully_genotyped = length(which((A_results_stats$concordance == TRUE | A_results_stats$concordance == FALSE)
                                                      & A_results_stats$caller == "GangSTR" )),
    A_num_eh_loci_fully_genotyped = length(which((A_results_stats$concordance == TRUE | A_results_stats$concordance == FALSE)
                                                 & A_results_stats$caller == "EH" )),
    A_num_concordant_loci = length(which(A_results_stats$concordance == TRUE )),
    A_num_hipstr_concordant_loci = length(which(A_results_stats$caller == "HipSTR" & A_results_stats$concordance == TRUE)),
    A_num_gangstr_concordant_loci = length(which(A_results_stats$caller == "GangSTR" & A_results_stats$concordance == TRUE)),
    A_num_eh_concordant_loci = length(which(A_results_stats$caller == "EH" & A_results_stats$concordance == TRUE)),
    A_num_discordant_loci = length(which(A_results_stats$concordance == FALSE )),
    A_num_hipstr_discordant_loci = length(which(A_results_stats$caller == "HipSTR" & A_results_stats$concordance == FALSE)),
    A_num_gangstr_discordant_loci = length(which(A_results_stats$caller == "GangSTR" & A_results_stats$concordance == FALSE)),
    A_num_eh_discordant_loci = length(which(A_results_stats$caller == "EH" & A_results_stats$concordance == FALSE)),
    A_num_loci_passing_filter_missing_gt = length(which(is.na(A_results_stats$concordance) & A_results_stats$caller != "no call")),
    A_num_hipstr_loci_passing_filter_missing_gt = length(which(A_results_stats$caller == "HipSTR" & is.na(A_results_stats$concordance))),
    A_num_gangstr_loci_passing_filter_missing_gt = length(which(A_results_stats$caller == "GangSTR" & is.na(A_results_stats$concordance))),
    A_num_eh_loci_passing_filter_missing_gt = length(which(A_results_stats$caller == "EH" & is.na(A_results_stats$concordance)))) %>%
  
  mutate(
    A_frac_num_passing_filters_over_A_num_loci_panel = A_num_loci_passing_filters / A_num_loci_panel,
    A_frac_loci_failing_filters_over_A_num_loci_panel = A_num_loci_failing_filters / A_num_loci_panel,
    A_frac_num_fully_genotyped_over_A_num_loci_panel = A_num_loci_fully_genotyped / A_num_loci_panel,
    A_frac_num_fully_genotyped_over_A_num_passing_filters = A_num_loci_fully_genotyped / A_num_loci_passing_filters,
    A_frac_num_caller_hipstr_over_A_num_passing_filters = A_num_loci_caller_hipstr / A_num_loci_passing_filters,
    A_frac_num_caller_gangstr_over_A_num_passing_filters = A_num_loci_caller_gangstr / A_num_loci_passing_filters,
    A_frac_num_caller_eh_over_A_num_passing_filters = A_num_loci_caller_eh / A_num_loci_passing_filters,
    A_frac_loci_passing_filter_missing_gt_over_A_num_passing_filters = A_num_loci_passing_filter_missing_gt / A_num_loci_passing_filters,
    A_frac_discordant_over_A_sum_concordant_discordant = A_num_discordant_loci / (A_num_concordant_loci + A_num_discordant_loci),
    A_frac_hipstr_discordant_over_A_sum_hipstr_concordant_discordant = 
      A_num_hipstr_discordant_loci / (A_num_hipstr_concordant_loci + A_num_hipstr_discordant_loci),
    A_frac_gangstr_discordant_over_A_sum_gangstr_concordant_discordant = 
      A_num_gangstr_discordant_loci / (A_num_gangstr_concordant_loci + A_num_gangstr_discordant_loci),
    A_frac_eh_discordant_over_A_sum_eh_concordant_discordant = 
      A_num_eh_discordant_loci / (A_num_eh_concordant_loci + A_num_eh_discordant_loci))


all_loci_stats <- cbind(nonA_stats, A_stats)

print(all_loci_stats)

# save concordance results in RDS
saveRDS(results, file = paste0(sub('\\.rds$', '', args[1]),"_",config$sample$trio,"_concordance.rds"))

# save concordance statistics in TSV file
write.table(all_loci_stats, file = paste0(sub('\\.rds$', '', args[1]),"_",config$sample$trio,"_concordance_stats.tsv"), quote = F, row.names = F)