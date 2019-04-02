#!/bin/bash

# Description: Use this script to assess the quality of cleaned (trimmed) reads. Sequence files of the forward and reverse reads must be indicated with "R1" and "R2" in the file names, respectively.
#
# Outputs: This script generates a report containing the lengths of all reads in each sequence file.
#
# Software requirements: 
#
# Usage: bash assess_cleaned_read_quality.sh [full path to the directory containing the demultiplexed/raw fastq/fastq.gz files] [name of directory containing cleaned reads]

#######################
# Specify directories #
#######################
rawReadsDir=$1

projDir=$(dirname $rawReadsDir)

resultsDir=${projDir}/processing_results
mkdir -p $resultsDir

cleanDirname=$2
cleanDir=${resultsDir}/${cleanDirname}

#####################################################################################################
# Create an array of the names that are used to identify the orientation of reads in a sequence file#
#####################################################################################################
orient=(R1 R2)

####################################################################################################################
# Create a single string of the patterns that must be removed from the file name to extract the name of the sample #
####################################################################################################################
fnameSubstr=${resultsDir}/filename_substrings_to_remove.txt

substr=($(cat ${fnameSubstr}))

############################################################
# Calculate lengths of the cleaned reads and write to a file
############################################################
for f in ${cleanDir}/*_cutadapt.fastq.gz
do
    fname=$(basename $f)
    
    for o in "${orient[@]}"
    do
        for str in "${substr[@]}"
        do
            pat=${str/_${o}_*fastq/_${o}_cutadapt.fastq}
            fname2=${fname/.gz/}
            sample=${fname2/${pat}/}
        done
    done
    
    gzcat $f | awk '{ if (NR%4==2) print length($1) }' > ${resultsDir}/${sample}_read_lengths.txt
done

########################################################
# Plot a histogram of the lengths of the cleaned reads #
########################################################
source activate /anaconda2

Rscript ~/Google_Drive/2018_ASU_Research/lizards/scripts/plot_hist_read_lengths.R $resultsDir

##############################################
# Assess the quality of the cleaning step(s) #
##############################################
Rscript ~/Google_Drive/2018_ASU_Research/lizards/scripts/calc_stats_plots_cleaned_reads.R $projDir mapped_reads_count.txt $cleanDirname

source deactivate /anaconda2