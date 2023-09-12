
#!/share/apps/r/4.3.1/gcc/bin/Rscript

# usage: Rscript join_str_calls_nextflow.r ${panel CSV} ${sampleID} ${sex} ${trio} ${hipstr} ${expansionhunter} ${gangstr} ${bedtools}

# Panel CSV is a csv that contains info directly from the shiny server for all loci in panel
# Panel CSV has 1-start, fully closed coordinates
# Requires panel CSV to be in the following order:
# row.number,seqnames,start,end,width,period.size,motif,motif.family

# HipSTR can't distinguish the sex of samples when performing joint calling
# This results in 2 alleles called on chromosomes X and Y in males when only one allele is present
# From the two alleles given in the H_GB field, we will choose the one associated with the higher H_MALLREADS count as the "true" allele
# If the true allele is in the allele 2 position initially in H_GB and H_MALLREADS, H_GB_1 and H_MALLREADS_1 will be swapped with H_GB_2 and H_MALLREADS_2
# This ensures that allele 1 in the final table is always the allele with most reads on the sex chromosomes
# We will replace H_REPDIFF_2 on sex chromosomes with NA to indicate allele 2 should not be considered for concordance
# If all samples are male, HipSTR will run with haploid calling for chrX and Y and this operation will not be necessary


library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = T)

# load list of panel loci and metadata
panel_loci <- as.data.frame(read.csv(args[1], header = F))
# add column names to panel loci data frame
colnames(panel_loci) <- c("name","chr","start","end","width","period.size","motif","motif.family")

# sample ID
sampleID <- args[2]

# sample sex (M or F)
sex <- args[3]

# trio ID
trio <- args[4]

#sample type (typeA or typeB)
sampleType <- args[5]

# load results table for each caller and bedtools
# results table for HipSTR from QUERY_HIPSTR
HipSTR <- read.delim(args[6])
# results table for Expansion Hunter from QUERY_EH
EH <- read.delim(args[7])
# results table for GangSTR from QUERY_GANGSTR
GangSTR <- read.delim(args[8])
# table of bedtools coverage from BEDTOOLS_COVERAGE
bedtools <- read.delim(args[9])


# make function for obtaining total read depth from MALLREADS, ALLREADS, and ENCLREADS
readdepth <- function(x) {sum(suppressWarnings(as.numeric(x))[seq(0, length(x), by=2)])}



# generate data frame containing all metadata and call information for the sample
results <- 
  # add sample ID and trio ID as columns to the panel loci data frame
  panel_loci %>% mutate(sampleID = sampleID, sex = sex, trio = trio, sampleType = sampleType) %>%
  # add column for motif length (since period size occasionally doesn't match)
  mutate(motifLength = nchar(motif.family)) %>%
  # left join caller and bedtools data files to panel loci data frame
  left_join(., HipSTR, by = "name") %>% left_join(., EH, by = "name") %>% 
  left_join(., GangSTR, by = c("chr","start")) %>% left_join(., bedtools, by = c("name")) %>%
  # replace any EH genotype that has B_depth of 0 with "./." on autosomes and "." on sex chromosomes
  mutate(E_GT = case_when(B_depth == 0 & sex == "M" & chr != "chrX" & chr != "chrY" ~ "./.",
                          B_depth == 0 & sex == "F" ~ "./.",
                          B_depth == 0 & sex == "M" & (chr == "chrX" | chr == "chrY") ~ ".", .default = E_GT)) %>%
  # replace any B_depth of 0 or B_depth with no EH call with NA for accurate calculation of mean depth
  mutate(B_depth = case_when(B_depth == 0 ~ NA, E_GT == "./." | E_GT == "." ~ NA, .default = B_depth)) %>%
  # separate fields giving base pair/repeat difference from reference into column for each allele
  separate(H_GB,c("H_GB_1_TEMP","H_GB_2_TEMP"), sep="\\|", fill="right", remove=F) %>%
  separate(E_REPCN, c("E_REPCN_1","E_REPCN_2"), sep="/", fill="right", remove=F) %>%
  separate(G_REPCN, c("G_REPCN_1","G_REPCN_2"), sep=",", fill="right", remove=F) %>%
  # split ADSP field to get spanning read depth for Expansion Hunter alleles
  separate(E_ADSP, c("E_ADSP_1","E_ADSP_2"), sep="/", fill="right", remove=F) %>%
  # convert columns from character to numeric where necessary
  mutate_at(c("H_Q","H_DP","H_DFLANKINDEL","H_DSTUTTER", "H_GB_1_TEMP", "H_GB_2_TEMP", "H_GLDIFF",
              "E_REPCN_1","E_REPCN_2","E_ADSP_1","E_ADSP_2",
              "G_Q","G_DP","G_REPCN_1","G_REPCN_2"), as.numeric) %>%
  # add columns to get repeat count difference from reference
  mutate(E_REPDIFF_1 = E_REPCN_1-E_REF, E_REPDIFF_2 = E_REPCN_2-E_REF) %>%
  mutate(G_REPDIFF_1 = G_REPCN_1-G_REF, G_REPDIFF_2 = G_REPCN_2-G_REF) %>%
  # add columns to get fraction of reads that have stutter or flank indels
  mutate(H_STUTTERFRAC = H_DSTUTTER/H_DP, H_FLINDELFRAC = H_DFLANKINDEL/H_DP)

# split MALLREADS and match to each allele from HipSTR call to get spanning read depth (maximum-likelihood-based)
results$H_MALLREADS_1_TEMP <- unlist(as.numeric(mapply(function(x,y){ sub(paste0("^",x,"\\|(.*)"), "\\1", y[grep(paste0("^",x,"\\|"),y)])}, x=results$H_GB_1_TEMP,
                                                       y=strsplit(results$H_MALLREADS,";"), SIMPLIFY=FALSE, USE.NAMES=FALSE)))
results$H_MALLREADS_2_TEMP <- unlist(as.numeric(mapply(function(x,y){ sub(paste0("^",x,"\\|(.*)"), "\\1", y[grep(paste0("^",x,"\\|"),y)])}, x=results$H_GB_2_TEMP,
                                                       y=strsplit(results$H_MALLREADS,";"), SIMPLIFY=FALSE, USE.NAMES=FALSE)))

# split ALLREADS and match to each allele from HipSTR call to get spanning read depth (alignment-based)
results$H_ALLREADS_1_TEMP <- unlist(as.numeric(mapply(function(x,y){ sub(paste0("^",x,"\\|(.*)"), "\\1", y[grep(paste0("^",x,"\\|"),y)])}, x=results$H_GB_1_TEMP,
                                                      y=strsplit(results$H_ALLREADS,";"), SIMPLIFY=FALSE, USE.NAMES=FALSE)))
results$H_ALLREADS_2_TEMP <- unlist(as.numeric(mapply(function(x,y){ sub(paste0("^",x,"\\|(.*)"), "\\1", y[grep(paste0("^",x,"\\|"),y)])}, x=results$H_GB_2_TEMP,
                                                      y=strsplit(results$H_ALLREADS,";"), SIMPLIFY=FALSE, USE.NAMES=FALSE)))

# Calculate repeat difference from reference for HipSTR based on sample sex
if (sex == "M") {
  
  results <- results %>%
    # fill in permanent H_MALLREADS columns
    # if on a sex chromosome, fill allele 1 with the allele with the maximum read count
    mutate(H_MALLREADS_1 = 
             if_else((chr == "chrX" | chr == "chrY") & H_MALLREADS_1_TEMP == pmax(H_MALLREADS_1_TEMP, H_MALLREADS_2_TEMP, na.rm = T), 
                     H_MALLREADS_1_TEMP, H_MALLREADS_2_TEMP),
           H_MALLREADS_2 =
             if_else((chr == "chrX" | chr == "chrY") & H_MALLREADS_1_TEMP == pmax(H_MALLREADS_1_TEMP, H_MALLREADS_2_TEMP, na.rm = T), 
                     H_MALLREADS_2_TEMP, H_MALLREADS_1_TEMP)) %>%
    # fill in permanent H_ALLREADS columns
    # if on a sex chromosome, fill allele 1 with the allele with the maximum read count
    mutate(H_ALLREADS_1 = 
             if_else((chr == "chrX" | chr == "chrY") & H_ALLREADS_1_TEMP == pmax(H_ALLREADS_1_TEMP, H_ALLREADS_2_TEMP, na.rm = T), 
                     H_ALLREADS_1_TEMP, H_ALLREADS_2_TEMP),
           H_ALLREADS_2 =
             if_else((chr == "chrX" | chr == "chrY") & H_ALLREADS_1_TEMP == pmax(H_ALLREADS_1_TEMP, H_ALLREADS_2_TEMP, na.rm = T), 
                     H_ALLREADS_2_TEMP, H_ALLREADS_1_TEMP)) %>%
    # fill in permanent H_GB columns
    # if on a sex chromosome, fill allele 1 with the allele with the maximum read count
    mutate(H_GB_1 = if_else((chr == "chrX" | chr == "chrY") & H_MALLREADS_1_TEMP == pmax(H_MALLREADS_1_TEMP, H_MALLREADS_2_TEMP, na.rm = T), 
                            H_GB_1_TEMP, H_GB_2_TEMP),
           H_GB_2 = if_else((chr == "chrX" | chr == "chrY") & H_MALLREADS_1_TEMP == pmax(H_MALLREADS_1_TEMP, H_MALLREADS_2_TEMP, na.rm = T), 
                            H_GB_2_TEMP, H_GB_1_TEMP)) %>%
    # calculate repeat difference
    # fill in NA for allele 2 on sex chromosomes to indicate it should not be considered a true call
    mutate(H_REPDIFF_1 = H_GB_1/motifLength,
           H_REPDIFF_2 = case_when((chr == "chrX" | chr == "chrY") ~ NA, .default = H_GB_2/motifLength)) %>%
    # calculate repeat count
    # fill in NA for allele 2 on sex chromosomes to indicate it should not be considered a true call
    mutate(H_REPCN_1 = (width + H_GB_1)/motifLength,
           H_REPCN_2 = case_when((chr == "chrX" | chr == "chrY") ~ NA, .default = (width + H_GB_2)/motifLength))
  
} else {
  
  # fill in all permanent columns with their corresponding temporary column
  # since sample is female, chrX has 2 alleles and can be filled as autosomes
  # fill in NA for any chrY calls to indicate it should not be considered a true call
  results <- results %>% 
    mutate(H_MALLREADS_1 = H_MALLREADS_1_TEMP, H_MALLREADS_2 = H_MALLREADS_2_TEMP,
           H_ALLREADS_1 = H_ALLREADS_1_TEMP, H_ALLREADS_2 = H_ALLREADS_2_TEMP,
           H_GB_1 = H_GB_1_TEMP, H_GB_2 = H_GB_2_TEMP) %>%
    mutate(H_REPDIFF_1 = case_when(chr == "chrY" ~ NA, .default = H_GB_1/motifLength), 
           H_REPDIFF_2 = case_when(chr == "chrY" ~ NA, .default = H_GB_2/motifLength)) %>%
    mutate(H_REPCN_1 = case_when(chr == "chrY" ~ NA, .default = (width + H_GB_1)/motifLength), 
           H_REPCN_2 = case_when(chr == "chrY" ~ NA, .default = (width + H_GB_2)/motifLength))
  
}

# remove temporary columns
results <- results %>% select(-c(H_GB_1_TEMP, H_GB_2_TEMP, H_MALLREADS_1_TEMP, H_MALLREADS_2_TEMP, H_ALLREADS_1_TEMP, H_ALLREADS_2_TEMP))


# split ENCLREADS and match to each allele from GangSTR call to get spanning read depth
results$G_ENCLREADS_1 <- unlist(as.numeric(mapply(function(x,y){ sub(paste0("^",x,",(.*)"), "\\1", y[grep(paste0("^",x,","),y)])}, x=results$G_REPCN_1,
                                                  y=strsplit(results$G_ENCLREADS,"\\|"), SIMPLIFY=FALSE, USE.NAMES=FALSE)))
results$G_ENCLREADS_2 <- unlist(as.numeric(mapply(function(x,y){ sub(paste0("^",x,",(.*)"), "\\1", y[grep(paste0("^",x,","),y)])}, x=results$G_REPCN_2,
                                                  y=strsplit(results$G_ENCLREADS,"\\|"), SIMPLIFY=FALSE, USE.NAMES=FALSE)))

# Calculate VAF
results <- results %>%
  # Find total spanning read depth in HipSTR by summing read depths in MALLREADS
  mutate(H_MALLREADS_SUM = sapply(strsplit(results$H_MALLREADS,";|\\|"), readdepth)) %>%
  mutate(H_ALLREADS_SUM = sapply(strsplit(results$H_ALLREADS,";|\\|"), readdepth)) %>%
  # Divide by individual read depth to get VAF
  mutate(H_MVAF_1 = H_MALLREADS_1/H_MALLREADS_SUM) %>% mutate(H_MVAF_2 = H_MALLREADS_2/H_MALLREADS_SUM) %>%
  mutate(H_AVAF_1 = H_ALLREADS_1/H_ALLREADS_SUM) %>% mutate(H_AVAF_2 = H_ALLREADS_2/H_ALLREADS_SUM) %>%
  mutate(E_VAF_1 = E_ADSP_1/B_depth) %>% mutate(E_VAF_2 = E_ADSP_2/B_depth) %>%
  # Find total spanning read depth in GangSTR by summing read depths in ENCLREADS
  mutate(G_ENCLREADS_SUM = sapply(strsplit(results$G_ENCLREADS,",|\\|"), readdepth)) %>%
  # Divide by individual read depth to get VAF
  mutate(G_VAF_1 = G_ENCLREADS_1/G_ENCLREADS_SUM) %>% mutate(G_VAF_2 = G_ENCLREADS_2/G_ENCLREADS_SUM)

# Save data frame in RDS file with sample name
saveRDS(results, file=paste0(sampleID,".rds"))