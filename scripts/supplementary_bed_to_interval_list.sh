#!/bin/bash

# SBATCH --time=100

# Usage: supp_bed_to_interval_list.sh [supplementary panel BED] [output file full path] [reference genome dictionary file]
# output file must end in .interval_list

# supplementary panel BED is a bed that contains the coordinates for the targets of the supplementary capture panel
# supplementary panel BED has 0-start, half-open coordinates
# requires supplementary panel BED to be tab-separated in the following order:
# [chr] [start] [end]
# Interval list coordinates are 0-start, half-open

module load jdk/11.0.9

java -jar /scratch/projects/evronylab/bin/gatk/picard.jar BedToIntervalList \
          I=$1 \
          O=$2 \
          SD=$3