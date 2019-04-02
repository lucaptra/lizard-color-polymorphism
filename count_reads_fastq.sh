#!/bin/bash

# Description: Use this script to count the number of reads that are in every fastq/fastq.gz file in a specified directory.
#
# Output: This script generates a tab-delimited file in which each line stores the name of the sequence file followed by the number of reads that are in that file.
#
# Software requirements: Standard BASH commands.
#
# Usage: bash count_reads_fastq.sh

#########################################
# Request values of variables from user #
#########################################
echo Please provide the full path to the directory containing the fastq/fastq.gz files for which you would like to count the number of reads: # ask user for path name of directory
read rawReadsDir

echo What would you like to call the output file? Please provide the file extension, too \(".txt" is recommended\). Don\'t include a directory path. # ask user for desired name of output file
read outFname

###########################
# Create output directory #
###########################
projDir=$(dirname $rawReadsDir) # get project directory

resultsDir=${projDir}/processing_results # set name of directory in which all outputs of pipeline will be written
mkdir -p $resultsDir

######################
# Create output file #
######################
out=${$resultsDir}/${outFname} # specify the full path of the output file

if [ -f "$out" ] # delete the output file from a previous run, if it exists
then
    rm $out
fi

echo -e filename'\t'#_reads > $out # add column names to the file

for file in ${rawReadsDir}/*R1*.fastq* # count the number of reads in each fastq/fastq.gz file in the directory
do
    if [ $file = *.gz ] # if the file is gzipped
    then
        gzcat $file | echo ${file}'\t'$((`wc -l` / 4)) >> $out
    else # if the file is uncompressed
        cat $file | echo $file'\t'$((`wc -l` / 4)) >> $out
    fi
done

echo "Number of reads tabulated for all fastq/fastq.gz files in your directory! Read counts written to file '"${out}"'". # print message when tasks are done

###############################################################
# Generate reports summarizing the number of reads per sample #
###############################################################
source activate /anaconda2

Rscript ~/Google_Drive/2018_ASU_Research/lizards/scripts/calc_stats_read_counts.R $out # write HTML table and plot histogram of number of fragments sequenced per sample

source deactivate /anaconda2