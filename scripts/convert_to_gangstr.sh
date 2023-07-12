
#!/bin/bash

# usage convert_to_gangstr.sh [Panel CSV] [full path to GangSTR input file]

# Panel CSV is a csv that contains info directly from the shiny server for all loci in panel
# Panel CSV has 1-start, fully closed coordinates
# requires panel CSV to be in the following order:
# row.number,seqnames,start,end,width,period.size,motif,motif.family

awk -F',' -e '{print $2"\t"$3"\t"$4"\t"$6"\t"$7}' $1 > $2

# Output tab-delimited format: [chr] [start] [end] [period size] [motif]
# Coordinates are 1-start, fully-closed