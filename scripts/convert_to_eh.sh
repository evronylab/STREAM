
#!/bin/bash

# Usage: convert_to_eh.sh [CSV of targeted microsatellites] [list of ExpansionHunter exclusion loci] [full path to eh_to_json.R] [full path to ExpansionHunter json output file]

# CSV of targeted microsatellites is a CSV that contains information about all loci in the target capture panel
# The CSV has 1-start, fully closed coordinates
# requires CSV to be in the following order:
# row.number,seqnames,start,end,width,period.size,motif,motif.family

# List of ExpansionHunter exclusion loci is a text file containing one locus name per line
# Locus name is in the format MS-#####
# The number in each name corresponds to the locus's row number in STRATIFY
# These loci cause ExpansionHunter to crash due to too many Ns in the flanks in the reference genome

# Temp file is in the following tab-delimited format:
# [LocusId] [LocusStructure] [ReferenceRegion] [VariantType]
# Make temp file in current directory so it can be accessed by R
# Acript subtracts 1 from shiny server start coordinate to make 0-based, half open (BED format)
# Final output is JSON file for ExpansionHunter 

TMPFILE=$(mktemp -p .)
OUT=$(mktemp -p .)

echo -e "LocusId\tLocusStructure\tReferenceRegion\tVariantType" > ${TMPFILE}
awk -F',' -e '{print $1"\t("$7")\*\t"$2":"($3-1)"-"$4"\tRepeat"}' $1 >> ${TMPFILE}

# Remove loci that cause ExpansionHunter to abort (no way to force program to skip these loci)
# Excluded loci cause error message "Error: too many Ns in nearby area"
grep -vwf $2 ${TMPFILE} > ${OUT}

# Send output to R script for conversion to JSON format
Rscript $3 ${OUT} $4

# Delete temp files when finished
rm ${TMPFILE}
rm ${OUT}