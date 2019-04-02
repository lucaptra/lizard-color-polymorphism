#!/bin/bash

# Description: Use this script to run FastQC and MultiQC in order to assess the quality of reads in every fastq/fastq.gz file in a specified directory. Sequence files containing the forward and reverse reads must be indicated with "R1" and "R2" in the file names, respectively.
#
# Outputs: This script generates one FastQC report for each sequence file and then summarizes the reports for each read orientation (i.e., forward or reverse). Documentation available for FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/) and MultiQC (http://multiqc.info/docs/).
#
# Software requirements: FastQC and MultiQC must be installed (e.g., "conda install -c bioconda fastqc -y"; "conda install -c bioconda multiqc -y").
#
# Usage: bash assess_read_quality_fastqc.sh

##########################################
# Request values of variables from user #
##########################################
echo Please provide the full path to the directory containing your fastq/fastq.gz files: # ask user for path name of directory
read rawReadsDir

echo Where would you like to write the reports from FastQC? Please provide a name (excluding path address) for the directory. # ask user for name of output directory
read fastqcDirname

###########################
# Create output directory #
###########################
projDir=$(dirname $rawReadsDir) # get project directory

resultsDir=${projDir}/processing_results # set name of directory in which all pipeline outputs will be written
fastqcDir=${resultsDir}/${fastqcDirname} # specify path name of FastQC output directory

for dir in ${fastqcDir} ${fastqcDir}/multiqc_reports_R1 ${fastqcDir}/multiqc_reports_R2 # create output directories to store reports from FastQC and MultiQC
do
	rm -Rf $dir # delete outputs from previous run, if they already exist
	mkdir -p $dir # create directory if it doesn't already exist
done

##############
# Run FastQC #
##############
for file in ${rawReadsDir}/*.fastq* # write quality control report for each fastq/fastq.gz file in the directory
do
    fastqc $file --outdir ${fastqcDir}
done

###############
# Run MultiQC #
###############
multiqc ${fastqcDir}/*.zip --ignore *_R2* --outdir ${fastqcDir}/multiqc_reports_R1 -f -p # summarize FastQC statistics for all forward reads in a single report
multiqc ${fastqcDir}/*.zip --ignore *_R1* --outdir ${fastqcDir}/multiqc_reports_R2 -f -p # summarize FastQC statistics for all reverse reads in a single report

echo FastQC reports are written for all fastq/fastq.gz files in your directory! Results are summarized across samples and written to ${fastqcDir}. # print message when tasks are done

#########################################################################
# Create HTML report summarizing quality control statistics from FastQC #
#########################################################################
source activate /anaconda2

Rscript ~/Google_Drive/2018_ASU_Research/lizards/scripts/summarize_fastqc_reports.R $fastqcDir # create HTML tables from MultiQC reports

source deactivate /anaconda2