#!/share/apps/r/4.3.1/gcc/bin/Rscript

# usage: Rscript usat_genotype_filtering.r [RDS output of usat_calling.nf] [YAML config]

# RDS is output of join_str_samples_nextflow.r with HipSTR, GangSTR, Expansion Hunter, and bedtools data for all loci and all samples in tibble format
# The RDS is the final output of the Nextflow pipeline
# config file is a YAML file with parameter settings for each caller and ensemble calling

library(dplyr)
library(tidyr)
library(stringr)
library(configr)

# load command line arguments
args <- commandArgs(trailingOnly = T)

# load results data frame from RDS
results <- readRDS(args[1])

# Load YAML config file
config <- read.config(file = args[2])

nonA_loci <- results %>%
  filter(motif.family != "A")

A_loci <- results %>%
  filter(motif.family == "A")


## filtering parameters

# HipSTR

# minimum HipSTR call quality (VCF field Q)
nonA.min.qual.hipstr <- config$HipSTR$nonA.min.qual
A.min.qual.hipstr <- config$HipSTR$A.min.qual
# minimum total number of reads used to make call in HipSTR (sum of read counts in VCF field MALLREADS)
nonA.min.mallreads.hipstr.typeA <- config$HipSTR$nonA.min.mallreads.typeA
nonA.min.mallreads.hipstr.typeB <- config$HipSTR$nonA.min.mallreads.typeB
A.min.mallreads.hipstr.typeA <- config$HipSTR$A.min.mallreads.typeA
A.min.mallreads.hipstr.typeB <- config$HipSTR$A.min.mallreads.typeB
# minimum number of reads supporting the called allele in HipSTR (read count for allele matching genotype call in VCF field MALLREADS)
nonA.min.allele.mallreads.hipstr.typeA <- config$HipSTR$nonA.min.allele.mallreads.typeA
nonA.min.allele.mallreads.hipstr.typeB <- config$HipSTR$nonA.min.allele.mallreads.typeB
A.min.allele.mallreads.hipstr.typeA <- config$HipSTR$A.min.allele.mallreads.typeA
A.min.allele.mallreads.hipstr.typeB <- config$HipSTR$A.min.allele.mallreads.typeB
# minimum variant allele fraction (VAF) for the called allele in HipSTR (MALLREADS for allele / sum of MALLREADS)
nonA.min.mvaf.hipstr <- config$HipSTR$nonA.min.mvaf
A.min.mvaf.hipstr <- config$HipSTR$A.min.mvaf
# minimum total number of reads used to make call in HipSTR (sum of read counts in VCF field ALLREADS)
nonA.min.allreads.hipstr.typeA <- config$HipSTR$nonA.min.allreads.typeA
nonA.min.allreads.hipstr.typeB <- config$HipSTR$nonA.min.allreads.typeB
A.min.allreads.hipstr.typeA <- config$HipSTR$A.min.allreads.typeA
A.min.allreads.hipstr.typeB <- config$HipSTR$A.min.allreads.typeB
# minimum number of reads supporting the called allele in HipSTR (read count for allele matching genotype call in VCF field allreads)
nonA.min.allele.allreads.hipstr.typeA <- config$HipSTR$nonA.min.allele.allreads.typeA
nonA.min.allele.allreads.hipstr.typeB <- config$HipSTR$nonA.min.allele.allreads.typeB
A.min.allele.allreads.hipstr.typeA <- config$HipSTR$A.min.allele.allreads.typeA
A.min.allele.allreads.hipstr.typeB <- config$HipSTR$A.min.allele.allreads.typeB
# minimum variant allele fraction (VAF) for the called allele in HipSTR (allreads for allele / sum of ALLREADS)
nonA.min.avaf.hipstr <- config$HipSTR$nonA.min.avaf
A.min.avaf.hipstr <- config$HipSTR$A.min.avaf
# maximum fraction of stutter reads
nonA.max.stutter.frac.hipstr <- config$HipSTR$nonA.max.stutter.frac
A.max.stutter.frac.hipstr <- config$HipSTR$A.max.stutter.frac
# maximum fraction of reads with flank indels
nonA.max.flankindel.frac.hipstr <- config$HipSTR$nonA.max.flankindel.frac
A.max.flankindel.frac.hipstr <- config$HipSTR$A.max.flankindel.frac
# minimum difference in likelihood between call and next-best call (GLDIFF)
nonA.min.gldiff.hipstr.typeA <- config$HipSTR$nonA.min.gldiff.typeA
nonA.min.gldiff.hipstr.typeB <- config$HipSTR$nonA.min.gldiff.typeB
A.min.gldiff.hipstr.typeA <- config$HipSTR$A.min.gldiff.typeA
A.min.gldiff.hipstr.typeB <- config$HipSTR$A.min.gldiff.typeB
# minimum average HipSTR call quality (VCF field Q) across all samples
nonA.min.mean.qual.hipstr <- config$HipSTR$nonA.min.mean.qual
A.min.mean.qual.hipstr <- config$HipSTR$A.min.mean.qual
# minimum average total number of reads (MALLREADS) used to make call in HipSTR across all samples
nonA.min.mean.mallreads.hipstr.typeA <- config$HipSTR$nonA.min.mean.mallreads.typeA
nonA.min.mean.mallreads.hipstr.typeB <- config$HipSTR$nonA.min.mean.mallreads.typeB
A.min.mean.mallreads.hipstr.typeA <- config$HipSTR$A.min.mean.mallreads.typeA
A.min.mean.mallreads.hipstr.typeB <- config$HipSTR$A.min.mean.mallreads.typeB
# minimum average total number of reads (ALLREADS) used to make call in HipSTR across all samples
nonA.min.mean.allreads.hipstr.typeA <- config$HipSTR$nonA.min.mean.allreads.typeA
nonA.min.mean.allreads.hipstr.typeB <- config$HipSTR$nonA.min.mean.allreads.typeB
A.min.mean.allreads.hipstr.typeA <- config$HipSTR$A.min.mean.allreads.typeA
A.min.mean.allreads.hipstr.typeB <- config$HipSTR$A.min.mean.allreads.typeB
# minimum fraction of total samples that have a call that passed all sample-level filters
nonA.min.sample.frac.hipstr <- config$HipSTR$nonA.min.sample.frac
A.min.sample.frac.hipstr <- config$HipSTR$A.min.sample.frac

# GangSTR

# minimum GangSTR call quality (VCF field Q)
nonA.min.qual.gangstr <- config$GangSTR$nonA.min.qual
A.min.qual.gangstr <- config$GangSTR$A.min.qual
# minimum total number of spanning reads used to make call in GangSTR (sum of read counts in VCF field ENCLREADS)
nonA.min.total.reads.gangstr.typeA <- config$GangSTR$nonA.min.total.reads.typeA
nonA.min.total.reads.gangstr.typeB <- config$GangSTR$nonA.min.total.reads.typeB
A.min.total.reads.gangstr.typeA <- config$GangSTR$A.min.total.reads.typeA
A.min.total.reads.gangstr.typeB <- config$GangSTR$A.min.total.reads.typeB
# minimum number of reads supporting the called allele in GangSTR (read count for allele matching genotype call in VCF field ENCLREADS)
nonA.min.allele.reads.gangstr.typeA <- config$GangSTR$nonA.min.allele.reads.typeA
nonA.min.allele.reads.gangstr.typeB <- config$GangSTR$nonA.min.allele.reads.typeB
A.min.allele.reads.gangstr.typeA <- config$GangSTR$A.min.allele.reads.typeA
A.min.allele.reads.gangstr.typeB <- config$GangSTR$A.min.allele.reads.typeB
# minimum variant allele fraction (VAF) for the called allele in GangSTR (ENCLREADS for allele / sum of ENCLREADS)
nonA.min.vaf.gangstr <- config$GangSTR$nonA.min.vaf
A.min.vaf.gangstr <- config$GangSTR$A.min.vaf
# minimum average GangSTR call quality (VCF field Q) across all samples
nonA.min.mean.qual.gangstr <- config$GangSTR$nonA.min.mean.qual
A.min.mean.qual.gangstr <- config$GangSTR$A.min.mean.qual
# minimum average total number of reads used to make call in GangSTR across all samples
nonA.min.mean.total.reads.gangstr.typeA <- config$GangSTR$nonA.min.mean.total.reads.typeA
nonA.min.mean.total.reads.gangstr.typeB <- config$GangSTR$nonA.min.mean.total.reads.typeB
A.min.mean.total.reads.gangstr.typeA <- config$GangSTR$A.min.mean.total.reads.typeA
A.min.mean.total.reads.gangstr.typeB <- config$GangSTR$A.min.mean.total.reads.typeB
# minimum fraction of total samples that have a call that passed all sample-level filters
nonA.min.sample.frac.gangstr <- config$GangSTR$nonA.min.sample.frac
A.min.sample.frac.gangstr <- config$GangSTR$A.min.sample.frac

# Expansion Hunter

# minimum number of spanning reads counted by bedtools coverage at each locus
nonA.min.total.reads.eh.typeA <- config$EH$nonA.min.total.reads.typeA
nonA.min.total.reads.eh.typeB <- config$EH$nonA.min.total.reads.typeB
A.min.total.reads.eh.typeA <- config$EH$A.min.total.reads.typeA
A.min.total.reads.eh.typeB <- config$EH$A.min.total.reads.typeB
# minimum number of reads supporting the called allele in EH (read count for allele matching genotype call in VCF field ADSP)
nonA.min.allele.reads.eh.typeA <- config$EH$nonA.min.allele.reads.typeA
nonA.min.allele.reads.eh.typeB <- config$EH$nonA.min.allele.reads.typeB
A.min.allele.reads.eh.typeA <- config$EH$A.min.allele.reads.typeA
A.min.allele.reads.eh.typeB <- config$EH$A.min.allele.reads.typeB
# minimum variant allele fraction (VAF) for the called allele in EH (ADSP for allele / bedtools coverage spanning reads)
nonA.min.vaf.eh <- config$EH$nonA.min.vaf
A.min.vaf.eh <- config$EH$A.min.vaf
# maximum variant allele fraction (VAF) for the called allele in EH (ADSP for allele / bedtools coverage spanning reads)
nonA.max.vaf.eh <- config$EH$nonA.max.vaf
A.max.vaf.eh <- config$EH$A.max.vaf
# minimum average locus coverage (LC) in EH
nonA.min.lc.eh <- config$EH$nonA.min.lc
A.min.lc.eh <- config$EH$A.min.lc
# minimum average number of spanning reads counted by bedtools coverage at each locus across all samples
nonA.min.mean.total.reads.eh.typeA <- config$EH$nonA.min.mean.total.reads.typeA
nonA.min.mean.total.reads.eh.typeB <- config$EH$nonA.min.mean.total.reads.typeB
A.min.mean.total.reads.eh.typeA <- config$EH$A.min.mean.total.reads.typeA
A.min.mean.total.reads.eh.typeB <- config$EH$A.min.mean.total.reads.typeB
# minimum fraction of total samples that have a call that passed all sample-level filters
nonA.min.sample.frac.eh <- config$EH$nonA.min.sample.frac
A.min.sample.frac.eh <- config$EH$A.min.sample.frac

# ensemble

# 1st preferred caller (HipSTR, GangSTR, or EH)
nonA.caller.one <- config$ensemble$nonA.caller.one
A.caller.one <- config$ensemble$A.caller.one
# 2nd preferred caller (HipSTR, GangSTR, or EH)
nonA.caller.two <- config$ensemble$nonA.caller.two
A.caller.two <- config$ensemble$A.caller.two
# 3rd preferred caller (HipSTR, GangSTR, or EH)
nonA.caller.three <- config$ensemble$nonA.caller.three
A.caller.three <- config$ensemble$A.caller.three


filtering_function <-
  function(subset_data,
           min.qual.hipstr, min.mallreads.hipstr.typeA, min.mallreads.hipstr.typeB, min.allreads.hipstr.typeA, min.allreads.hipstr.typeB,
           min.allele.mallreads.hipstr.typeA, min.allele.mallreads.hipstr.typeB, min.allele.allreads.hipstr.typeA, min.allele.allreads.hipstr.typeB, min.mvaf.hipstr, min.avaf.hipstr,
           max.stutter.frac.hipstr, max.flankindel.frac.hipstr, min.gldiff.hipstr.typeA, min.gldiff.hipstr.typeB,
           min.mean.qual.hipstr, min.mean.mallreads.hipstr.typeA, min.mean.mallreads.hipstr.typeB, min.mean.allreads.hipstr.typeA, min.mean.allreads.hipstr.typeB, min.sample.frac.hipstr,
           min.qual.gangstr, min.total.reads.gangstr.typeA, min.total.reads.gangstr.typeB, min.allele.reads.gangstr.typeA, min.allele.reads.gangstr.typeB, min.vaf.gangstr,
           min.mean.qual.gangstr, min.mean.total.reads.gangstr.typeA, min.mean.total.reads.gangstr.typeB, min.sample.frac.gangstr,
           min.total.reads.eh.typeA, min.total.reads.eh.typeB, min.allele.reads.eh.typeA, min.allele.reads.eh.typeB, min.vaf.eh, max.vaf.eh, min.lc.eh,
           min.mean.total.reads.eh.typeA, min.mean.total.reads.eh.typeB, min.sample.frac.eh,
           caller.one, caller.two, caller.three) {
    
    ## Sample-level filters
    
    # add column for each caller
    # fill with TRUE if sample passes all filters for that caller and chromosome
    # fill with FALSE if sample fails any filter for that caller
    
    if((caller.one == "no_caller")) {
      stop("caller.one cannot be set to 'no_caller'")
    }
    
    if((caller.two == "no_caller" & caller.three != "no_caller")){
      stop("If caller.two == 'no_caller', caller.three must also be set to 'no_caller'")
    }
    
    subset_data <- subset_data %>%
      # HipSTR
      mutate(HipSTR_passsample =
               # if on chrX or chrY and male
               case_when(sex == "M" & (chr == "chrX" | chr == "chrY") &
                           # if there is a genotype call
                           H_GT != "." & !is.na(H_GT) &
                           # quality
                           H_Q >= min.qual.hipstr & !is.na(H_Q) &
                           # total read depth
                           if_else(sampleType == "typeA", H_MALLREADS_SUM >= min.mallreads.hipstr.typeA, H_MALLREADS_SUM >= min.mallreads.hipstr.typeB) & !is.na(H_MALLREADS_SUM) &
                           if_else(sampleType == "typeA", H_ALLREADS_SUM >= min.allreads.hipstr.typeA, H_ALLREADS_SUM >= min.allreads.hipstr.typeB) & !is.na(H_ALLREADS_SUM) &                       
                           # allele 1 read depth
                           if_else(sampleType == "typeA", H_MALLREADS_1 >= min.allele.mallreads.hipstr.typeA, H_MALLREADS_1 >= min.allele.mallreads.hipstr.typeB) & !is.na(H_MALLREADS_1) &
                           if_else(sampleType == "typeA", H_ALLREADS_1 >= min.allele.allreads.hipstr.typeA, H_ALLREADS_1 >= min.allele.allreads.hipstr.typeB) & !is.na(H_ALLREADS_1) &                       
                           # allele 1 VAF
                           H_MVAF_1 >= min.mvaf.hipstr & !is.na(H_MVAF_1) &
                           H_AVAF_1 >= min.avaf.hipstr & !is.na(H_AVAF_1) &                       
                           # stutter fraction
                           H_STUTTERFRAC <= max.stutter.frac.hipstr & !is.na(H_STUTTERFRAC) &
                           # flank indel fraction
                           H_FLINDELFRAC <= max.flankindel.frac.hipstr & !is.na(H_FLINDELFRAC) &
                           # likelihood difference
                           # GLDIFF is sometimes not reported for loci with genotypes
                           # We do not want to throw out those loci                       
                           (if_else(sampleType == "typeA", H_GLDIFF >= min.gldiff.hipstr.typeA, H_GLDIFF >= min.gldiff.hipstr.typeB) | is.na(H_GLDIFF)) ~ TRUE,
                         
                         # if on an autosome and male or on any chromosome other than chrY and female
                         ((sex == "M" & chr != "chrX" & chr != "chrY") | (sex == "F" & chr != "chrY")) &
                           # if there is a genotype call
                           H_GT != "." & !is.na(H_GT) &
                           # quality
                           H_Q >= min.qual.hipstr & !is.na(H_Q) &
                           # total read depth
                           if_else(sampleType == "typeA", H_MALLREADS_SUM >= min.mallreads.hipstr.typeA, H_MALLREADS_SUM >= min.mallreads.hipstr.typeB) & !is.na(H_MALLREADS_SUM) &
                           if_else(sampleType == "typeA", H_ALLREADS_SUM >= min.allreads.hipstr.typeA, H_ALLREADS_SUM >= min.allreads.hipstr.typeB) & !is.na(H_ALLREADS_SUM) &                       
                           # allele 1 read depth
                           if_else(sampleType == "typeA", H_MALLREADS_1 >= min.allele.mallreads.hipstr.typeA, H_MALLREADS_1 >= min.allele.mallreads.hipstr.typeB) & !is.na(H_MALLREADS_1) &
                           if_else(sampleType == "typeA", H_ALLREADS_1 >= min.allele.allreads.hipstr.typeA, H_ALLREADS_1 >= min.allele.allreads.hipstr.typeB) & !is.na(H_ALLREADS_1) &                       
                           # allele 2 read depth
                           if_else(sampleType == "typeA", H_MALLREADS_2 >= min.allele.mallreads.hipstr.typeA, H_MALLREADS_2 >= min.allele.mallreads.hipstr.typeB) & !is.na(H_MALLREADS_2) &
                           if_else(sampleType == "typeA", H_ALLREADS_2 >= min.allele.allreads.hipstr.typeA, H_ALLREADS_2 >= min.allele.allreads.hipstr.typeB) & !is.na(H_ALLREADS_2) &                       
                           # allele 1 VAF
                           H_MVAF_1 >= min.mvaf.hipstr & !is.na(H_MVAF_1) &
                           H_AVAF_1 >= min.avaf.hipstr & !is.na(H_AVAF_1) &                       
                           # allele 2 VAF
                           H_MVAF_2 >= min.mvaf.hipstr & !is.na(H_MVAF_2) &
                           H_AVAF_2 >= min.avaf.hipstr & !is.na(H_AVAF_2) &                       
                           # stutter fraction
                           H_STUTTERFRAC <= max.stutter.frac.hipstr & !is.na(H_STUTTERFRAC) &
                           # flank indel fraction
                           H_FLINDELFRAC <= max.flankindel.frac.hipstr & !is.na(H_FLINDELFRAC) &
                           # likelihood difference
                           # GLDIFF is sometimes not reported for loci with genotypes
                           # we do not want to throw out those loci
                           # different thresholds based on sampleType
                           (if_else(sampleType == "typeA", H_GLDIFF >= min.gldiff.hipstr.typeA, H_GLDIFF >= min.gldiff.hipstr.typeB) | is.na(H_GLDIFF)) ~ TRUE,
                         
                         # if any filter failed
                         .default = FALSE)) %>%
      
      # GangSTR
      mutate(GangSTR_passsample =
               # if on chrX or chr Y and male
               case_when(sex == "M" & (chr == "chrX" | chr == "chrY") &
                           # if there is a genotype call
                           G_GT != "." & !is.na(G_GT) &
                           # quality
                           G_Q >= min.qual.gangstr & !is.na(G_Q) &
                           # total read depth
                           if_else(sampleType == "typeA", G_ENCLREADS_SUM >= min.total.reads.gangstr.typeA, G_ENCLREADS_SUM >= min.total.reads.gangstr.typeB) & !is.na(G_ENCLREADS_SUM) &
                           # allele 1 read depth
                           if_else(sampleType == "typeA", G_ENCLREADS_1 >= min.allele.reads.gangstr.typeA, G_ENCLREADS_1 >= min.allele.reads.gangstr.typeB) & !is.na(G_ENCLREADS_1) &
                           # allele 1 VAF
                           G_VAF_1 >= min.vaf.gangstr & !is.na(G_VAF_1) ~ TRUE,
                         
                         # if on an autosome and male or on any chromosome other than chrY and female
                         ((sex == "M" & chr != "chrX" & chr != "chrY") | (sex == "F" & chr != "chrY")) &
                           # if there is a genotype call
                           G_GT != "." & !is.na(G_GT) &
                           # quality
                           G_Q >= min.qual.gangstr & !is.na(G_Q) &
                           # total read depth
                           if_else(sampleType == "typeA", G_ENCLREADS_SUM >= min.total.reads.gangstr.typeA, G_ENCLREADS_SUM >= min.total.reads.gangstr.typeB) & !is.na(G_ENCLREADS_SUM) &
                           # allele 1 read depth
                           if_else(sampleType == "typeA", G_ENCLREADS_1 >= min.allele.reads.gangstr.typeA, G_ENCLREADS_1 >= min.allele.reads.gangstr.typeB) & !is.na(G_ENCLREADS_1) &
                           # allele 2 read depth
                           if_else(sampleType == "typeA", G_ENCLREADS_2 >= min.allele.reads.gangstr.typeA, G_ENCLREADS_2 >= min.allele.reads.gangstr.typeB) & !is.na(G_ENCLREADS_2) &
                           # allele 1 VAF
                           G_VAF_1 >= min.vaf.gangstr & !is.na(G_VAF_1) &
                           # allele 2 VAF
                           G_VAF_2 >= min.vaf.gangstr & !is.na(G_VAF_2) ~ TRUE,
                         
                         # if any filter failed
                         .default = FALSE)) %>%
      
      # Expansion Hunter
      mutate(EH_passsample =
               # if on chrX or chrY and male
               case_when(sex == "M" & (chr == "chrX" | chr == "chrY") &
                           # if there is a genotype call
                           E_GT != "." & !is.na(E_GT) &
                           # alleles supported by spanning reads only
                           (E_SO == "SPANNING" | E_SO == "SPANNING/SPANNING") & !is.na(E_SO) &
                           # total read depth (bedtools)
                           if_else(sampleType == "typeA", B_depth >= min.total.reads.eh.typeA, B_depth >= min.total.reads.eh.typeB) & !is.na(B_depth) &
                           # allele 1 read depth
                           if_else(sampleType == "typeA", E_ADSP_1 >= min.allele.reads.eh.typeA, E_ADSP_1 >= min.allele.reads.eh.typeB) & !is.na(E_ADSP_1) &
                           # allele 1 VAF
                           E_VAF_1 >= min.vaf.eh & E_VAF_1 <= max.vaf.eh & !is.na(E_VAF_1) &
                           # average locus coverage
                           E_LC >= min.lc.eh & !is.na(E_LC) ~ TRUE,
                         
                         # if on an autosome and male or on any chromosome other than chrY and female
                         ((sex == "M" & chr != "chrX" & chr != "chrY") | (sex == "F" & chr != "chrY")) &
                           # if there is a genotype call
                           E_GT != "./." & !is.na(E_GT) &
                           # alleles supported by spanning reads only
                           (E_SO == "SPANNING" | E_SO == "SPANNING/SPANNING") & !is.na(E_SO) &
                           # total read depth (bedtools)
                           if_else(sampleType == "typeA", B_depth >= min.total.reads.eh.typeA, B_depth >= min.total.reads.eh.typeB) & !is.na(B_depth) &
                           # allele 1 read depth
                           if_else(sampleType == "typeA", E_ADSP_1 >= min.allele.reads.eh.typeA, E_ADSP_1 >= min.allele.reads.eh.typeB) & !is.na(E_ADSP_1) &
                           # allele 2 read depth
                           if_else(sampleType == "typeA", E_ADSP_2 >= min.allele.reads.eh.typeA, E_ADSP_2 >= min.allele.reads.eh.typeB) & !is.na(E_ADSP_2) &
                           # allele 1 VAF
                           E_VAF_1 >= min.vaf.eh & E_VAF_1 <= max.vaf.eh & !is.na(E_VAF_1) &
                           # allele 2 VAF
                           E_VAF_2 >= min.vaf.eh & E_VAF_2 <= max.vaf.eh & !is.na(E_VAF_2) &
                           # average locus coverage
                           E_LC >= min.lc.eh & !is.na(E_LC) ~ TRUE,
                         
                         # if any filter failed
                         .default = FALSE))
    
    ## Locus-level filters
    
    # get number of samples from table
    num_samples = n_distinct(subset_data$sampleID)
    
    # Calculate and apply locus-level filters
    qual_avg <- subset_data %>% group_by(name) %>%
      summarise(# average quality
        mean_qual_H = mean(H_Q, na.rm = T),
        mean_qual_G = mean(G_Q, na.rm = T))
    depth_avg <- subset_data %>% group_by(name, sampleType) %>%
      summarise(# average total read depth
        mean_mdepth_H = mean(H_MALLREADS_SUM, na.rm = T),
        mean_adepth_H = mean(H_ALLREADS_SUM, na.rm = T),
        mean_depth_G = mean(G_ENCLREADS_SUM, na.rm = T),
        mean_depth_E = mean(B_depth, na.rm = T))
    depth_qual_avg <- left_join(depth_avg, qual_avg, by = "name") %>%
      mutate(
        #HipSTR filter
        locus_means_H = case_when(sampleType == "typeA" & mean_qual_H >= min.mean.qual.hipstr &
                                    mean_mdepth_H >= min.mean.mallreads.hipstr.typeA &
                                    mean_adepth_H >= min.mean.allreads.hipstr.typeA ~ T,
                                  sampleType == "typeB" & mean_qual_H >= min.mean.qual.hipstr &
                                    mean_mdepth_H >= min.mean.mallreads.hipstr.typeB &
                                    mean_adepth_H >= min.mean.allreads.hipstr.typeB ~ T,
                                  .default = F),
        # GangSTR filter
        locus_means_G = case_when(sampleType == "typeA" & mean_qual_G >= min.mean.qual.gangstr & mean_depth_G >= min.mean.total.reads.gangstr.typeA ~ T,
                                  sampleType == "typeB" & mean_qual_G >= min.mean.qual.gangstr & mean_depth_G >= min.mean.total.reads.gangstr.typeB ~ T,
                                  .default = F),
        # Expansion Hunter filter
        locus_means_E = case_when(sampleType == "typeA" & mean_depth_E >= min.mean.total.reads.eh.typeA ~ T
                                  sampleType == "typeB" & mean_depth_E >= min.mean.total.reads.eh.typeB ~ T
                                  .default = F)) %>%
      # replace NA values in locus means columns with FALSE
      mutate(locus_means_H = if_else(is.na(locus_means_H) == T, F, locus_means_H),
             locus_means_G = if_else(is.na(locus_means_G) == T, F, locus_means_G),
             locus_means_E = if_else(is.na(locus_means_E) == T, F, locus_means_E))
    
    
    # Calculate the fraction of samples in the data frame that passed the sample-level filters for each locus and caller
    # If the sample passed sample-level filters, the column [HipSTR|GangSTR|EH]_passsample in the results data frame == TRUE
    # The number of samples passing sample-level filters at each locus is found using summarise(sum([HipSTR|GangSTR|EH]_passsample == TRUE))
    # Then that number is divided by the total number of samples to get the fraction of samples passing at each locus
    # If the resulting fraction is > min.sample.frac.[hipstr|gangstr|eh], sample_frac_filter_[H|G|E] == TRUE for that locus
    # If the resulting fraction is < min.sample.frac.[hipstr|gangstr|eh], sample_frac_filter_[H|G|E] == FALSE for that locus
    sample_frac_filter <- subset_data %>% select(c(name, sampleID, HipSTR_passsample, GangSTR_passsample, EH_passsample)) %>%
      group_by(name) %>%
      # HipSTR
      summarise(sample_frac_filter_H = if_else(sum(HipSTR_passsample)/num_samples >= min.sample.frac.hipstr, T, F),
                # GangSTR
                sample_frac_filter_G = if_else(sum(GangSTR_passsample)/num_samples >= min.sample.frac.gangstr, T, F),
                # EH
                sample_frac_filter_E = if_else(sum(EH_passsample)/num_samples >= min.sample.frac.eh, T, F))
    
    
    # Join to main table
    subset_data <- left_join(subset_data, sample_frac_filter, by = "name") %>% left_join(., depth_qual_avg, by = c("name", "sampleType"))
    
    # Use Boolean logic to decide whether to keep locus
    subset_data <- subset_data %>% group_by(name) %>%
      mutate(HipSTR_passlocus = all(locus_means_H) & sample_frac_filter_H,
             GangSTR_passlocus = all(locus_means_G) & sample_frac_filter_G,
             EH_passlocus = all(locus_means_E) & sample_frac_filter_E)
    
    # Label caller.one calls, and set remaining to "no call"
    subset_data <- subset_data %>% mutate(caller = if_else(!!as.name(paste0(caller.one,"_passlocus")) == T,caller.one,"no call"))
    
    # Label caller.two calls
    if(caller.two != "no_caller"){
      subset_data <- subset_data %>% mutate(caller=if_else(caller == "no call" & !!as.name(paste0(caller.two,"_passlocus")) == T,caller.two,caller))
    }
    
    # Label caller.three calls
    if(caller.three != "no_caller"){
      subset_data <- subset_data %>% mutate(caller=if_else(caller == "no call" & !!as.name(paste0(caller.three,"_passlocus")) == T,caller.three,caller))
    }
    
    
    # Extract genotypes for all loci with a passing call
    
    ## REPEAT DIFFERENCE
    
    # Fill in allele_1
    # Genotypes for caller.one calls, and set remaining to NA
    subset_data <- subset_data %>% mutate(allele_1 = if_else(caller == caller.one & !!as.name(paste0(caller.one,"_passsample")) == T, !!as.name(paste0(substr(caller.one,1,1),"_REPDIFF_1")), NA_real_))
    # Genotypes for caller.two calls
    if(caller.two != "no_caller") {
      subset_data <- subset_data %>% mutate(allele_1 = if_else(caller == caller.two & !!as.name(paste0(caller.two,"_passsample")) == T, !!as.name(paste0(substr(caller.two,1,1),"_REPDIFF_1")), allele_1))
    }
    # Genotypes for caller.three calls
    if(caller.three != "no_caller") {
      subset_data <- subset_data %>% mutate(allele_1 = if_else(caller == caller.three & !!as.name(paste0(caller.three,"_passsample")) == T, !!as.name(paste0(substr(caller.three,1,1),"_REPDIFF_1")), allele_1))
    }
    
    # Fill in allele_2
    # Genotypes for caller.one calls, and set remaining to NA
    subset_data <- subset_data %>% mutate(allele_2 = if_else(caller == caller.one & !!as.name(paste0(caller.one,"_passsample")) == T, !!as.name(paste0(substr(caller.one,1,1),"_REPDIFF_2")), NA_real_))
    # Genotypes for caller.two calls
    if(caller.two != "no_caller") {
      subset_data <- subset_data %>% mutate(allele_2 = if_else(caller == caller.two & !!as.name(paste0(caller.two,"_passsample")) == T, !!as.name(paste0(substr(caller.two,1,1),"_REPDIFF_2")), allele_2))
    }
    # Genotypes for caller.three calls
    if(caller.three != "no_caller") {
      subset_data <- subset_data %>% mutate(allele_2 = if_else(caller == caller.three & !!as.name(paste0(caller.three,"_passsample")) == T, !!as.name(paste0(substr(caller.three,1,1),"_REPDIFF_2")), allele_2))
    }
    
    ## REPEAT COUNT
    # Fill in allele_1
    # Genotypes for caller.one calls, and set remaining to NA
    subset_data <- subset_data %>% mutate(allele_1_repcn = if_else(caller == caller.one & !!as.name(paste0(caller.one,"_passsample")) == T, !!as.name(paste0(substr(caller.one,1,1),"_REPCN_1")), NA_real_))
    # Genotypes for caller.two calls
    if(caller.two != "no_caller") {
      subset_data <- subset_data %>% mutate(allele_1_repcn = if_else(caller == caller.two & !!as.name(paste0(caller.two,"_passsample")) == T, !!as.name(paste0(substr(caller.two,1,1),"_REPCN_1")), allele_1_repcn))
    }
    # Genotypes for caller.three calls
    if(caller.three != "no_caller") {
      subset_data <- subset_data %>% mutate(allele_1_repcn = if_else(caller == caller.three & !!as.name(paste0(caller.three,"_passsample")) == T, !!as.name(paste0(substr(caller.three,1,1),"_REPCN_1")), allele_1_repcn))
    }
    
    # Fill in allele_2
    # Genotypes for caller.one calls, and set remaining to NA
    subset_data <- subset_data %>% mutate(allele_2_repcn = if_else(caller == caller.one & !!as.name(paste0(caller.one,"_passsample")) == T, !!as.name(paste0(substr(caller.one,1,1),"_REPCN_2")), NA_real_))
    # Genotypes for caller.two calls
    if(caller.two != "no_caller") {
      subset_data <- subset_data %>% mutate(allele_2_repcn = if_else(caller == caller.two & !!as.name(paste0(caller.two,"_passsample")) == T, !!as.name(paste0(substr(caller.two,1,1),"_REPCN_2")), allele_2_repcn))
    }
    # Genotypes for caller.three calls
    if(caller.three != "no_caller") {
      subset_data <- subset_data %>% mutate(allele_2_repcn = if_else(caller == caller.three & !!as.name(paste0(caller.three,"_passsample")) == T, !!as.name(paste0(substr(caller.three,1,1),"_REPCN_2")), allele_2_repcn))
    }
    
    
    return(subset_data)
    
  } # END OF FILTERING FUNCTION

# Run filtering on non-poly-A loci
nonA_loci <-
  filtering_function(nonA_loci,
                     nonA.min.qual.hipstr, nonA.min.mallreads.hipstr.typeA, nonA.min.mallreads.hipstr.typeB, nonA.min.allreads.hipstr.typeA, nonA.min.allreads.hipstr.typeB,
                     nonA.min.allele.mallreads.hipstr.typeA, nonA.min.allele.mallreads.hipstr.typeB, nonA.min.allele.allreads.hipstr.typeA, nonA.min.allele.allreads.hipstr.typeB,
                     nonA.min.mvaf.hipstr, nonA.min.avaf.hipstr, 
                     nonA.max.stutter.frac.hipstr, nonA.max.flankindel.frac.hipstr, nonA.min.gldiff.hipstr.typeA, nonA.min.gldiff.hipstr.typeB,
                     nonA.min.mean.qual.hipstr, nonA.min.mean.mallreads.hipstr.typeA, nonA.min.mean.mallreads.hipstr.typeB, nonA.min.mean.allreads.hipstr.typeA, nonA.min.mean.allreads.hipstr.typeB, nonA.min.sample.frac.hipstr,
                     nonA.min.qual.gangstr, nonA.min.total.reads.gangstr.typeA, nonA.min.total.reads.gangstr.typeB, nonA.min.allele.reads.gangstr.typeA, nonA.min.allele.reads.gangstr.typeB,
                     nonA.min.vaf.gangstr, nonA.min.mean.qual.gangstr, nonA.min.mean.total.reads.gangstr.typeA, nonA.min.mean.total.reads.gangstr.typeB, nonA.min.sample.frac.gangstr,
                     nonA.min.total.reads.eh.typeA, nonA.min.total.reads.eh.typeB, nonA.min.allele.reads.eh.typeA, nonA.min.allele.reads.eh.typeB,
                     nonA.min.vaf.eh, nonA.max.vaf.eh, nonA.min.lc.eh,
                     nonA.min.mean.total.reads.eh.typeA, nonA.min.mean.total.reads.eh.typeB, nonA.min.sample.frac.eh,
                     nonA.caller.one, nonA.caller.two, nonA.caller.three)


# Run filtering on poly-A loci
A_loci <-
  filtering_function(A_loci,
                     A.min.qual.hipstr, A.min.mallreads.hipstr.typeA, A.min.mallreads.hipstr.typeB, A.min.allreads.hipstr.typeA, A.min.allreads.hipstr.typeB,
                     A.min.allele.mallreads.hipstr.typeA, A.min.allele.mallreads.hipstr.typeB, A.min.allele.allreads.hipstr.typeA, A.min.allele.allreads.hipstr.typeB,
                     A.min.mvaf.hipstr, A.min.avaf.hipstr, 
                     A.max.stutter.frac.hipstr, A.max.flankindel.frac.hipstr, A.min.gldiff.hipstr.typeA, A.min.gldiff.hipstr.typeB,
                     A.min.mean.qual.hipstr, A.min.mean.mallreads.hipstr.typeA, A.min.mean.mallreads.hipstr.typeB, A.min.mean.allreads.hipstr.typeA, A.min.mean.allreads.hipstr.typeB, A.min.sample.frac.hipstr,
                     A.min.qual.gangstr, A.min.total.reads.gangstr.typeA, A.min.total.reads.gangstr.typeB, A.min.allele.reads.gangstr.typeA, A.min.allele.reads.gangstr.typeB,
                     A.min.vaf.gangstr, A.min.mean.qual.gangstr, A.min.mean.total.reads.gangstr.typeA, A.min.mean.total.reads.gangstr.typeB, A.min.sample.frac.gangstr,
                     A.min.total.reads.eh.typeA, A.min.total.reads.eh.typeB, A.min.allele.reads.eh.typeA, A.min.allele.reads.eh.typeB,
                     A.min.vaf.eh, A.max.vaf.eh, A.min.lc.eh,
                     A.min.mean.total.reads.eh.typeA, A.min.mean.total.reads.eh.typeB, A.min.sample.frac.eh,
                     A.caller.one, A.caller.two, A.caller.three)


# Report which filter caused a locus to fail

error_report_function <-
  function(subset_data,
           min.qual.hipstr, min.mallreads.hipstr.typeA, min.mallreads.hipstr.typeB, min.allreads.hipstr.typeA, min.allreads.hipstr.typeB,
           min.allele.mallreads.hipstr.typeA, min.allele.mallreads.hipstr.typeB, min.allele.allreads.hipstr.typeA, min.allele.allreads.hipstr.typeB, min.mvaf.hipstr, min.avaf.hipstr,
           max.stutter.frac.hipstr, max.flankindel.frac.hipstr, min.gldiff.hipstr.typeA, min.gldiff.hipstr.typeB,
           min.mean.qual.hipstr, min.mean.mallreads.hipstr.typeA, min.mean.mallreads.hipstr.typeB, min.mean.allreads.hipstr.typeA, min.mean.allreads.hipstr.typeB, min.sample.frac.hipstr,
           min.qual.gangstr, min.total.reads.gangstr.typeA, min.total.reads.gangstr.typeB, min.allele.reads.gangstr.typeA, min.allele.reads.gangstr.typeB, min.vaf.gangstr,
           min.mean.qual.gangstr, min.mean.total.reads.gangstr.typeA, min.mean.total.reads.gangstr.typeB, min.sample.frac.gangstr,
           min.total.reads.eh.typeA, min.total.reads.eh.typeB, min.allele.reads.eh.typeA, min.allele.reads.eh.typeB, min.vaf.eh, max.vaf.eh, min.lc.eh,
           min.mean.total.reads.eh.typeA, min.mean.total.reads.eh.typeB, min.sample.frac.eh) {
    
    ## Sample-level filters
    
    # fill with name of filter if failed
    # fill with blank if passed
    # allele 2 filters are applied only to autosomes and female sex chromsomes
    
    subset_data <- subset_data %>%
      
      mutate(sample_failed_filters = "") %>%
      mutate(sample_failed_filters = if_else(H_GT != "." & !is.na(H_GT),
                                             sample_failed_filters,str_c(sample_failed_filters,"no.gt.hipstr",";"))) %>%
      mutate(sample_failed_filters = if_else(H_Q >= min.qual.hipstr & !is.na(H_Q),
                                             sample_failed_filters,str_c(sample_failed_filters,"min.qual.hipstr",";"))) %>%
      mutate(sample_failed_filters = case_when(sampleType == "typeA" & H_MALLREADS_SUM >= min.mallreads.hipstr.typeA & !is.na(H_MALLREADS_SUM) ~ sample_failed_filters,
                                               sampleType == "typeB" & H_MALLREADS_SUM >= min.mallreads.hipstr.typeB & !is.na(H_MALLREADS_SUM) ~ sample_failed_filters,
                                               .default = str_c(sample_failed_filters,"min.mallreads.hipstr",";"))) %>%
      mutate(sample_failed_filters = case_when(sampleType == "typeA" & H_ALLREADS_SUM >= min.allreads.hipstr.typeA & !is.na(H_ALLREADS_SUM) ~ sample_failed_filters,
                                               sampleType == "typeB" & H_ALLREADS_SUM >= min.allreads.hipstr.typeB & !is.na(H_ALLREADS_SUM) ~ sample_failed_filters,
                                               .default = str_c(sample_failed_filters,"min.allreads.hipstr",";"))) %>%
      mutate(sample_failed_filters = case_when(sampleType == "typeA" & H_MALLREADS_1 >= min.allele.mallreads.hipstr.typeA & !is.na(H_MALLREADS_1) ~ sample_failed_filters,
                                               sampleType == "typeB" & H_MALLREADS_1 >= min.allele.mallreads.hipstr.typeB & !is.na(H_MALLREADS_1) ~ sample_failed_filters,
                                               .default = str_c(sample_failed_filters,"min.allele.mallreads.hipstr_1",";"))) %>%
      mutate(sample_failed_filters = case_when(((sex == "M" & chr != "chrX" & chr != "chrY") | (sex == "F" & chr != "chrY")) & sampleType == "typeA" &
                                               H_MALLREADS_2 >= min.allele.mallreads.hipstr.typeA & !is.na(H_MALLREADS_2) ~ sample_failed_filters,
                                               ((sex == "M" & chr != "chrX" & chr != "chrY") | (sex == "F" & chr != "chrY")) & sampleType == "typeB" &
                                                 H_MALLREADS_2 >= min.allele.mallreads.hipstr.typeB & !is.na(H_MALLREADS_2) ~ sample_failed_filters,
                                               .default = str_c(sample_failed_filters,"min.allele.mallreads.hipstr_2",";"))) %>%
      mutate(sample_failed_filters = case_when(sampleType == "typeA" & H_ALLREADS_1 >= min.allele.allreads.hipstr.typeA & !is.na(H_ALLREADS_1) ~ sample_failed_filters,
                                               sampleType == "typeB" & H_ALLREADS_1 >= min.allele.allreads.hipstr.typeB & !is.na(H_ALLREADS_1) ~ sample_failed_filters,
                                               .default = str_c(sample_failed_filters,"min.allele.allreads.hipstr_1",";"))) %>%
      mutate(sample_failed_filters = case_when(((sex == "M" & chr != "chrX" & chr != "chrY") | (sex == "F" & chr != "chrY")) & sampleType == "typeA" &
                                                 H_ALLREADS_2 >= min.allele.allreads.hipstr.typeA & !is.na(H_ALLREADS_2) ~ sample_failed_filters,
                                               ((sex == "M" & chr != "chrX" & chr != "chrY") | (sex == "F" & chr != "chrY")) & sampleType == "typeB" &
                                                 H_ALLREADS_2 >= min.allele.allreads.hipstr.typeB & !is.na(H_ALLREADS_2) ~ sample_failed_filters,
                                               .default = str_c(sample_failed_filters,"min.allele.allreads.hipstr_2",";"))) %>%
      mutate(sample_failed_filters = if_else(H_MVAF_1 >= min.mvaf.hipstr & !is.na(H_MVAF_1),
                                             sample_failed_filters,str_c(sample_failed_filters,"min.mvaf.hipstr_1",";"))) %>%
      mutate(sample_failed_filters = if_else(((sex == "M" & chr != "chrX" & chr != "chrY") | (sex == "F" & chr != "chrY")) &
                                               H_MVAF_2 >= min.mvaf.hipstr & !is.na(H_MVAF_2),
                                             sample_failed_filters,str_c(sample_failed_filters,"min.mvaf.hipstr_2",";"))) %>%
      mutate(sample_failed_filters = if_else(H_AVAF_1 >= min.avaf.hipstr & !is.na(H_AVAF_1),
                                             sample_failed_filters,str_c(sample_failed_filters,"min.avaf.hipstr_1",";"))) %>%
      mutate(sample_failed_filters = if_else(((sex == "M" & chr != "chrX" & chr != "chrY") | (sex == "F" & chr != "chrY")) &
                                               H_AVAF_2 >= min.avaf.hipstr & !is.na(H_AVAF_2),
                                             sample_failed_filters,str_c(sample_failed_filters,"min.avaf.hipstr_2",";"))) %>%
      mutate(sample_failed_filters = if_else(H_STUTTERFRAC <= max.stutter.frac.hipstr & !is.na(H_STUTTERFRAC),
                                             sample_failed_filters,
                                             str_c(sample_failed_filters,"max.stutter.frac.hipstr",";"))) %>%
      mutate(sample_failed_filters = if_else(H_FLINDELFRAC <= max.flankindel.frac.hipstr & !is.na(H_FLINDELFRAC),
                                             sample_failed_filters,
                                             str_c(sample_failed_filters,"max.flankindel.frac.hipstr",";"))) %>%
      mutate(sample_failed_filters = case_when(sampleType == "typeA" & (H_GLDIFF >= min.gldiff.hipstr.typeA | is.na(H_GLDIFF)) ~ sample_failed_filters,
                                               sampleType == "typeB" & (H_GLDIFF >= min.gldiff.hipstr.typeB | is.na(H_GLDIFF)) ~ sample_failed_filters,
                                               .default = str_c(sample_failed_filters,"min.gldiff.hipstr",";"))) %>%
      mutate(sample_failed_filters = if_else(G_GT != "." & !is.na(G_GT),
                                             sample_failed_filters,str_c(sample_failed_filters,"no.gt.gangstr",";"))) %>%
      mutate(sample_failed_filters = if_else(G_Q >= min.qual.gangstr & !is.na(G_Q),
                                             sample_failed_filters,str_c(sample_failed_filters,"min.qual.gangstr",";"))) %>%
      mutate(sample_failed_filters = case_when(sampleType == "typeA" & G_ENCLREADS_SUM >= min.total.reads.gangstr.typeA & !is.na(G_ENCLREADS_SUM) ~ sample_failed_filters,
                                               sampleType == "typeB" & G_ENCLREADS_SUM >= min.total.reads.gangstr.typeB & !is.na(G_ENCLREADS_SUM) ~ sample_failed_filters,
                                               .default = str_c(sample_failed_filters,"min.total.reads.gangstr",";"))) %>%
      mutate(sample_failed_filters = if_else(sampleType == "typeA" & G_ENCLREADS_1 >= min.allele.reads.gangstr.typeA & !is.na(G_ENCLREADS_1) ~ sample_failed_filters,
                                             sampleType == "typeB" & G_ENCLREADS_1 >= min.allele.reads.gangstr.typeB & !is.na(G_ENCLREADS_1) ~ sample_failed_filters,
                                             .default = str_c(sample_failed_filters,"min.allele.reads.gangstr_1",";"))) %>%
      mutate(sample_failed_filters = if_else(((sex == "M" & chr != "chrX" & chr != "chrY") | (sex == "F" & chr != "chrY")) & sampleType == "typeA" &
                                               G_ENCLREADS_2 >= min.allele.reads.gangstr.typeA & !is.na(G_ENCLREADS_2)  ~ sample_failed_filters,
                                             ((sex == "M" & chr != "chrX" & chr != "chrY") | (sex == "F" & chr != "chrY")) & sampleType == "typeB" &
                                               G_ENCLREADS_2 >= min.allele.reads.gangstr.typeB & !is.na(G_ENCLREADS_2)  ~ sample_failed_filters,
                                             .default = str_c(sample_failed_filters,"min.allele.reads.gangstr_2",";"))) %>%
      mutate(sample_failed_filters = if_else(G_VAF_1 >= min.vaf.gangstr & !is.na(G_VAF_1),
                                             sample_failed_filters,str_c(sample_failed_filters,"min.vaf.gangstr_1",";"))) %>%
      mutate(sample_failed_filters = if_else(((sex == "M" & chr != "chrX" & chr != "chrY") | (sex == "F" & chr != "chrY")) &
                                               G_VAF_2 >= min.vaf.gangstr & !is.na(G_VAF_2),
                                             sample_failed_filters,str_c(sample_failed_filters,"min.vaf.gangstr_2",";"))) %>%
      mutate(sample_failed_filters = if_else(E_GT != "." & !is.na(E_GT),
                                             sample_failed_filters,str_c(sample_failed_filters,"no.gt.eh",";"))) %>%
      mutate(sample_failed_filters = if_else(sampleType == "typeA" & B_depth >= min.total.reads.eh.typeA & !is.na(B_depth) ~ sample_failed_filters,
                                             sampleType == "typeB" & B_depth >= min.total.reads.eh.typeB & !is.na(B_depth) ~ sample_failed_filters,
                                             .default = str_c(sample_failed_filters,"min.total.reads.eh",";"))) %>%
      mutate(sample_failed_filters = if_else(sampleType == "typeA" & E_ADSP_1 >= min.allele.reads.eh.typeA & !is.na(E_ADSP_1) ~ sample_failed_filters,
                                             sampleType == "typeB" & E_ADSP_1 >= min.allele.reads.eh.typeB & !is.na(E_ADSP_1) ~ sample_failed_filters,
                                             .default = str_c(sample_failed_filters,"min.allele.reads.eh_1",";"))) %>%
      mutate(sample_failed_filters = if_else(((sex == "M" & chr != "chrX" & chr != "chrY") | (sex == "F" & chr != "chrY")) & sampleType == "typeA" &
                                               E_ADSP_2 >= min.allele.reads.eh.typeA & !is.na(E_ADSP_2) ~ sample_failed_filters,
                                             ((sex == "M" & chr != "chrX" & chr != "chrY") | (sex == "F" & chr != "chrY")) & sampleType == "typeB" &
                                               E_ADSP_2 >= min.allele.reads.eh.typeB & !is.na(E_ADSP_2) ~ sample_failed_filters,
                                             .default = str_c(sample_failed_filters,"min.allele.reads.eh_2",";"))) %>%
      mutate(sample_failed_filters = if_else(E_VAF_1 >= min.vaf.eh & E_VAF_1 <= max.vaf.eh & !is.na(E_VAF_1),
                                             sample_failed_filters,str_c(sample_failed_filters,"min.vaf.eh_1",";"))) %>%
      mutate(sample_failed_filters = if_else(((sex == "M" & chr != "chrX" & chr != "chrY") | (sex == "F" & chr != "chrY")) &
                                               E_VAF_2 >= min.vaf.eh & E_VAF_2 <= max.vaf.eh & !is.na(E_VAF_2),
                                             sample_failed_filters,str_c(sample_failed_filters,"min.vaf.eh_2",";"))) %>%
      mutate(sample_failed_filters = if_else(E_LC >= min.lc.eh & !is.na(E_LC),
                                             sample_failed_filters,str_c(sample_failed_filters,"min.lc.eh",";"))) %>%
      
      ## Locus-level filters
      
      # fill with name of filter if failed
      # fill with blank if passed
      
      mutate(locus_failed_filters = "") %>%
      mutate(locus_failed_filters = if_else(mean_qual_H >= min.mean.qual.hipstr & !is.na(mean_qual_H),
                                            locus_failed_filters,str_c(locus_failed_filters,"min.mean.qual.hipstr",";"))) %>%
      mutate(locus_failed_filters = if_else(mean_qual_G >= min.mean.qual.gangstr & !is.na(mean_qual_G),
                                            locus_failed_filters,str_c(locus_failed_filters,"min.mean.qual.gangstr",";"))) %>%
      mutate(locus_failed_filters = if_else(sample_frac_filter_H == T,
                                            locus_failed_filters,str_c(locus_failed_filters,"min.sample.frac.hipstr",";"))) %>%
      mutate(locus_failed_filters = if_else(sample_frac_filter_G == T,
                                            locus_failed_filters,str_c(locus_failed_filters,"min.sample.frac.gangstr",";"))) %>%
      mutate(locus_failed_filters = if_else(sample_frac_filter_E == T,
                                            locus_failed_filters,str_c(locus_failed_filters,"min.sample.frac.eh",";")))
  
    # table for checking mean read depth filters
    # handles typeA and typeB failures
    mean_depth_failed_filters <- subset_data %>% 
      mutate(mean_mdepth_H_pass_A = case_when(sampleType == "typeA" & mean_mdepth_H >= min.mean.mallreads.hipstr.typeA ~ T,
                                              sampleType == "typeB" ~ NA, .default = F),
             mean_mdepth_H_pass_B = case_when(sampleType == "typeB" & mean_mdepth_H >= min.mean.mallreads.hipstr.typeB ~ T,
                                              sampleType == "typeA" ~ NA, .default = F),
             mean_adepth_H_pass_A = case_when(sampleType == "typeA" & mean_adepth_H >= min.mean.allreads.hipstr.typeA ~ T,
                                              sampleType == "typeB" ~ NA, .default = F),
             mean_adepth_H_pass_B = case_when(sampleType == "typeB" & mean_adepth_H >= min.mean.allreads.hipstr.typeB ~ T,
                                              sampleType == "typeA" ~ NA, .default = F),
             mean_depth_G_pass_A = case_when(sampleType == "typeA" & mean_depth_G >= min.mean.total.reads.gangstr.typeA ~ T,
                                             sampleType == "typeB" ~ NA, .default = F),
             mean_depth_G_pass_B = case_when(sampleType == "typeB" & mean_depth_G >= min.mean.total.reads.gangstr.typeB ~ T,
                                             sampleType == "typeA" ~ NA, .default = F),
             mean_depth_E_pass_A = case_when(sampleType == "typeA" & mean_depth_E >= min.mean.total.reads.eh.typeA ~ T,
                                             sampleType == "typeB" ~ NA, .default = F),
             mean_depth_E_pass_B = case_when(sampleType == "typeB" & mean_depth_E >= min.mean.total.reads.eh.typeB ~ T,
                                             sampleType == "typeA" ~ NA, .default = F)) %>%
      group_by(name) %>% summarise(
        all_mean_mdepth_H_pass_A = all(mean_mdepth_H_pass_A == T, na.rm = T),
        all_mean_mdepth_H_pass_B = all(mean_mdepth_H_pass_B == T, na.rm = T),
        all_mean_adepth_H_pass_A = all(mean_adepth_H_pass_A == T, na.rm = T),
        all_mean_adepth_H_pass_B = all(mean_adepth_H_pass_B == T, na.rm = T),
        all_mean_depth_G_pass_A = all(mean_depth_G_pass_A == T, na.rm = T),
        all_mean_depth_G_pass_B = all(mean_depth_G_pass_B == T, na.rm = T),
        all_mean_depth_E_pass_A = all(mean_depth_E_pass_A == T, na.rm = T),
        all_mean_depth_E_pass_B = all(mean_depth_E_pass_B == T, na.rm = T)
      )
    
    # join to main table and fill in error column
    # remove extra columns when done
    subset_data <- left_join(subset_data, mean_depth_failed_filters, by = "name") %>%
      mutate(locus_failed_filters = if_else(all_mean_mdepth_H_pass_A == T, locus_failed_filters,str_c(locus_failed_filters,"min.mean.mallreads.hipstr.typeA",";"))) %>%
      mutate(locus_failed_filters = if_else(all_mean_mdepth_H_pass_B == T, locus_failed_filters,str_c(locus_failed_filters,"min.mean.mallreads.hipstr.typeB",";"))) %>%
      mutate(locus_failed_filters = if_else(all_mean_adepth_H_pass_A == T, locus_failed_filters,str_c(locus_failed_filters,"min.mean.allreads.hipstr.typeA",";"))) %>%
      mutate(locus_failed_filters = if_else(all_mean_adepth_H_pass_B == T, locus_failed_filters,str_c(locus_failed_filters,"min.mean.allreads.hipstr.typeB",";"))) %>%
      mutate(locus_failed_filters = if_else(all_mean_depth_G_pass_A == T, locus_failed_filters,str_c(locus_failed_filters,"min.mean.total.reads.gangstr.typeA",";"))) %>%
      mutate(locus_failed_filters = if_else(all_mean_depth_G_pass_B == T, locus_failed_filters,str_c(locus_failed_filters,"min.mean.total.reads.gangstr.typeB",";"))) %>%
      mutate(locus_failed_filters = if_else(all_mean_depth_E_pass_A == T, locus_failed_filters,str_c(locus_failed_filters,"min.mean.total.reads.eh.typeA",";"))) %>%
      mutate(locus_failed_filters = if_else(all_mean_depth_E_pass_B == T, locus_failed_filters,str_c(locus_failed_filters,"min.mean.total.reads.eh.typeB",";"))) %>%
      
      select(-c(all_mean_mdepth_H_pass_A,all_mean_mdepth_H_pass_B,all_mean_adepth_H_pass_A,all_mean_adepth_H_pass_B,
                all_mean_depth_G_pass_A,all_mean_depth_G_pass_B,all_mean_depth_E_pass_A,all_mean_depth_E_pass_B))
    

    return(subset_data)
    
  } # END OF ERROR FUNCTION


# Run error reporting on non-poly-A loci
nonA_loci <-
  error_report_function(nonA_loci,
                        nonA.min.qual.hipstr, nonA.min.mallreads.hipstr.typeA, nonA.min.mallreads.hipstr.typeB, nonA.min.allreads.hipstr.typeA, nonA.min.allreads.hipstr.typeB,
                        nonA.min.allele.mallreads.hipstr.typeA, nonA.min.allele.mallreads.hipstr.typeB, nonA.min.allele.allreads.hipstr.typeA, nonA.min.allele.allreads.hipstr.typeB,
                        nonA.min.mvaf.hipstr, nonA.min.avaf.hipstr, 
                        nonA.max.stutter.frac.hipstr, nonA.max.flankindel.frac.hipstr, nonA.min.gldiff.hipstr.typeA, nonA.min.gldiff.hipstr.typeB,
                        nonA.min.mean.qual.hipstr, nonA.min.mean.mallreads.hipstr.typeA, nonA.min.mean.mallreads.hipstr.typeB, nonA.min.mean.allreads.hipstr.typeA, nonA.min.mean.allreads.hipstr.typeB, nonA.min.sample.frac.hipstr,
                        nonA.min.qual.gangstr, nonA.min.total.reads.gangstr.typeA, nonA.min.total.reads.gangstr.typeB, nonA.min.allele.reads.gangstr.typeA, nonA.min.allele.reads.gangstr.typeB,
                        nonA.min.vaf.gangstr, nonA.min.mean.qual.gangstr, nonA.min.mean.total.reads.gangstr.typeA, nonA.min.mean.total.reads.gangstr.typeB, nonA.min.sample.frac.gangstr,
                        nonA.min.total.reads.eh.typeA, nonA.min.total.reads.eh.typeB, nonA.min.allele.reads.eh.typeA, nonA.min.allele.reads.eh.typeB,
                        nonA.min.vaf.eh, nonA.max.vaf.eh, nonA.min.lc.eh,
                        nonA.min.mean.total.reads.eh.typeA, nonA.min.mean.total.reads.eh.typeB, nonA.min.sample.frac.eh)


# Run error reporting on poly-A loci
A_loci <-
  error_report_function(A_loci,
                        A.min.qual.hipstr, A.min.mallreads.hipstr.typeA, A.min.mallreads.hipstr.typeB, A.min.allreads.hipstr.typeA, A.min.allreads.hipstr.typeB,
                        A.min.allele.mallreads.hipstr.typeA, A.min.allele.mallreads.hipstr.typeB, A.min.allele.allreads.hipstr.typeA, A.min.allele.allreads.hipstr.typeB,
                        A.min.mvaf.hipstr, A.min.avaf.hipstr, 
                        A.max.stutter.frac.hipstr, A.max.flankindel.frac.hipstr, A.min.gldiff.hipstr.typeA, A.min.gldiff.hipstr.typeB,
                        A.min.mean.qual.hipstr, A.min.mean.mallreads.hipstr.typeA, A.min.mean.mallreads.hipstr.typeB, A.min.mean.allreads.hipstr.typeA, A.min.mean.allreads.hipstr.typeB, A.min.sample.frac.hipstr,
                        A.min.qual.gangstr, A.min.total.reads.gangstr.typeA, A.min.total.reads.gangstr.typeB, A.min.allele.reads.gangstr.typeA, A.min.allele.reads.gangstr.typeB,
                        A.min.vaf.gangstr, A.min.mean.qual.gangstr, A.min.mean.total.reads.gangstr.typeA, A.min.mean.total.reads.gangstr.typeB, A.min.sample.frac.gangstr,
                        A.min.total.reads.eh.typeA, A.min.total.reads.eh.typeB, A.min.allele.reads.eh.typeA, A.min.allele.reads.eh.typeB,
                        A.min.vaf.eh, A.max.vaf.eh, A.min.lc.eh,
                        A.min.mean.total.reads.eh.typeA, A.min.mean.total.reads.eh.typeB, A.min.sample.frac.eh)


# Define chromosome sort order
chrOrder <-c(paste0("chr",1:22),"chrX", "chrY")

# Join non-A and poly-A loci and sort by chromosome and start coordinate
filtered_results <- rbind(nonA_loci, A_loci) %>%
  arrange(sampleID, factor(chr, chrOrder), start)

# Save genotypes list in RDS
saveRDS(filtered_results, file = paste0(sub('\\.rds$', '', args[1]),"_genotypes.rds"))
