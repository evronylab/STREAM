
#!/share/apps/r/4.2.0/gcc/bin/Rscript

# usage: Rscript eh_to_json.R [temp file from convert_to_eh.sh script] [full path to output json]
# Script will run automatically from convert_to_eh.sh
# Input file from convert_to_eh.sh is in the following tab-delimited format:
# [LocusId] [LocusStructure] [ReferenceRegion] [VariantType]

library(dplyr)
library(jsonlite)

# Import arguments from command line (in this case, from convert_to_EH.sh)
args <- commandArgs(trailingOnly = T)

# Read bed file of loci created in convert_to_EH.sh
eh_loci <- read.delim(args[1])

# Convert bed file to JSON format
# Use "pretty = T" to get proper formatting with new lines and indents
panel_json <- toJSON(eh_loci, pretty = T)

# Write JSON file to file name given to convert_to_EH.sh
write(panel_json, args[2])