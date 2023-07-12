#!/bin/bash

# SBATCH --time=100

# Usage: bed_to_interval_list.sh [Panel CSV] [Probe coordinates file] [Reference genome dictionary file]

# Panel CSV is a CSV that contains info directly from the shiny server for all loci in panel
# Panel CSV has 1-start, fully closed coordinates
# requires panel CSV to be in the following order:
# row.number,seqnames,start,end,width,period.size,motif,motif.family

# Probe coordinates file is a txt file that contains the coordinates of all regions covered by probes
# Probe coordinates file has 1-start, fully closed coordinates
# requires probe coordinates file to be tab-separated in the following order:
# [chr] [start] [end] [locus name]


module load jdk/11.0.9

USAT=$(mktemp -p .)
BAIT=$(mktemp -p .)

awk -F',' -e '{print $2"\t"($3-1)"\t"$4}' $1 > ${USAT}
awk -F'\t' -e '{print $1"\t"($2-1)"\t"$3}' $2 > ${BAIT}

# Make interval list of microsatellite coordinates
java -jar /scratch/projects/evronylab/bin/gatk/picard.jar BedToIntervalList \
          I=${USAT} \
          O=$(dirname $2)/$(basename $1 .csv).interval_list \
          SD=$3

# Make interval list of probe coordinates
java -jar /scratch/projects/evronylab/bin/gatk/picard.jar BedToIntervalList \
          I=${BAIT} \
          O=$(dirname $2)/$(basename $2 .txt).interval_list \
          SD=$3

rm ${USAT}
rm ${BAIT}

# Interval list coordinates are 0-start, half-open