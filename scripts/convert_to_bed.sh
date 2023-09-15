
#!/bin/bash

# usage convert_to_bed.sh [CSV of targeted microsatellites] [full path to bedtools input file]

# CSV of targeted microsatellites is a CSV that contains information about all loci in the target capture panel
# The CSV has 1-start, fully closed coordinates
# requires CSV to be in the following order:
# row.number,seqnames,start,end,width,period.size,motif,motif.family

awk -F',' -e '{print $2"\t"($3-1)"\t"$4"\t"$1}' $1 > $2

# Output tab-delimited BED file: [chr] [start] [end] [locus name]
# Coordinates are 0-start, half-open