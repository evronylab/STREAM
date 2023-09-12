#!/bin/bash

# SBATCH --time=100

# Usage: convert_to_interval_list.sh [Targets/probes TSV] [output file basename] [Reference genome dictionary file]

# Targets/probes TSV is a tab-separated file that contains the coordinates of the targeted loci or capture probes
# The TSV file should have 0-start, fully closed coordinates
# The TSV file must be tab-separated and have the following fields:
# [chr] [start coordinate] [end coordinate]

module load jdk/11.0.9

# Make interval list of coordinates for the targets or probes
java -jar /scratch/projects/evronylab/bin/gatk/picard.jar BedToIntervalList \
          I=$1 \
          O=$2.interval_list \
          SD=$3

# Interval list coordinates are 0-start, half-open