
#!/share/apps/r/4.3.1/gcc/bin/Rscript

#uusage: Rscript join_str_samples_nextflow.r ${RDS list} ${params.outputbasename}

library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = T)

# Initialize results list
results_list <- list()

# Load each RDS file and join data frame to list
for (x in head(args,-1)) {
  results_list[[x]] <- readRDS(x)
}

# Bind each data frame row-wise
# Final row count should equal number of samples * number of loci in panel
results <- as_tibble(do.call(rbind, results_list))

# Save list in RDS file
saveRDS(results, file = paste0(tail(args, 1),".rds"))