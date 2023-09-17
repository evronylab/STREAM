
#!/bin/bash

# usage convert_to_bed.sh [main microsatellites metadata CSV] [full path to bedtools input file]

# The microsatellites metadata CSV contains information about all loci to profile
# The CSV has 1-start, fully closed coordinates
# requires CSV to be in the following order:
# row.number,seqnames,start,end,width,period.size,motif,motif.family

awk -F',' -e '{print $2"\t"($3-1)"\t"$4"\t"$1}' $1 > $2

# Output tab-delimited BED file: [chr] [start] [end] [locus name]
# Coordinates are 0-start, half-open