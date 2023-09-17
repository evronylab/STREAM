
#!/bin/bash

# usage convert_to_hipstr.sh [main microsatellites metadata CSV] [full path to output HipSTR input file]
# Formats list into HipSTR region file

# The microsatellites metadata CSV contains information about all loci to profile
# The CSV has 1-start, fully closed coordinates
# requires CSV to be in the following order:
# row.number,seqnames,start,end,width,period.size,motif,motif.family

awk -F',' -e '{print $2"\t"$3"\t"$4"\t"$6"\t"$5/$6"\t"$1}' $1 > $2

# Output tab-delimited format: [chr] [start] [end] [period size] [number of repeats] [name]
# Coordinates are 1-start, fully-closed