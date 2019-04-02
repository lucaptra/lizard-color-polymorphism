#!/bin/bash

# Description: Use this script to calculate statistics for evaluating the quality of mapping to the reference genome.
#
# Outputs: This script generates a report containing various statistics for each SAM input file.
#
# Software requirements: 
#
# Usage: bash calc_mapping_stats_ngm.sh

#################
# Set variables #
#################
echo Provide a short, descriptive name for your data set. E.g., control-subset-0.1-s0.00-i0.65
read dataSet

echo What portion of the raw data was mapped for each sample?
read subsetPer

sample='N/A'
numReadsMerged='N/A'
mappedChr='N/A'
unmappedChr='N/A'
mapFracPerChr='N/A'
uniqReadsChr='N/A'
multiReadsChr='N/A'
mapqAvgChr='N/A'
numBases5XcovChr='N/A'
per5XcovChr='N/A'

##############################################
# Define parameter values to calculate stats #
##############################################
echo What was the sensitivity parameter \(-s\) used in the NextGenMap analysis? 
read sVal

echo What was the minimum identity parameter \(-i\) used in the NextGenMap analysis? 
read minIdent

#######################
# Specify directories #
#######################
echo Where are the outputs of NGM stored? Provide the full path name of the directory.
read ngmDir

projDir=$(dirname $ngmDir)
refDir=${projDir}/assemblies

######################################################################
# Check that the NGM analysis used the Anolis carolinensis reference #
######################################################################
echo Were the reads mapped to the Anolis carolinensis reference genome? yes or no.
read ref

if [ "$ref" == 'no' ]
then
    echo Sorry\, you cannot use this script for another reference genome.
    exit 0
fi

##########################################################################################################################
# Create literal string of chromosome names as query term to extract from sample BAM files reads that map to chromosomes #
##########################################################################################################################
chrNamesSearch='^'$(cat ${refDir}/chrNames.txt | sed -e 's/ /|^/g')

####################################
# BED files for specifying regions #
#####################################
bedChr=${refDir}/GCA_000090745.2_AnoCar2.0_genomic_chr.bed
bedScaf=${refDir}/GCA_000090745.2_AnoCar2.0_genomic_scaf.bed

################################################################
# Determine lengths of different regions of reference sequence #
################################################################
## Total sequenced length of reference
refSeq=${refDir}/Anolis_carolinensis.AnoCar2.0.dna_sm.toplevel.fa

refLen=$(awk '{s+=$0} END {printf "%d", s}' ${refDir}/Anolis_carolinensis.AnoCar2.0.dna_sm.toplevel_contig_lengths.txt)

###########################################################
## Calculate sequenced length of chromosomes in reference #
###########################################################
refLenChr=$(awk '{s+=$0} END {printf "%d", s}' ${refDir}/Anolis_carolinensis.AnoCar2.0.dna_sm.toplevel_contig_lengths_chr.txt)

#########################################################
## Calculate sequenced length of scaffolds in reference #
#########################################################
refLenScaf=$(awk '{s+=$0} END {printf "%d", s}' ${refDir}/Anolis_carolinensis.AnoCar2.0.dna_sm.toplevel_contig_lengths_scaf.txt)

################################################################
# Create output file to store mapping stats for set of samples #
################################################################
outTable=${ngmDir}/${dataSet/ /_}_mapping_stats_table.txt

rm -f $outTable
rm -f ${ngmDir}/*.bam.bai

###################################
# Set column names of output file #
###################################
header='data_set,sample,subset,s_value,min_identity_value,num_raw_reads,num_reads_merged,num_reads_mapped_chr,num_reads_mapped_scaf,num_reads_mapped_total,num_reads_unmapped_chr,num_reads_unmapped_scaf,%_reads_mapped_chr,%_reads_mapped_scaf,%_reads_mapped_total,num_mapped_reads_unique_chr,num_mapped_reads_unique_scaf,num_mapped_reads_multi_chr,num_mapped_reads_multi_scaf,mapped_reads_mapq_avg_chr,mapped_reads_mapq_avg_scaf,num_bases_>=5X_cov_chr,num_bases_>=5X_cov_scaf,num_bases_>=5X_cov_total,%_breadth_coverage_chr,%_breadth_coverage_scaf,%_breadth_coverage_AnoCar2.0'

###################################
# Add column names to output file #
###################################
echo ${header//,/$'\t'} >> $outTable

############################
# Ask for name of SAM file #
############################
echo What is the file name of the SAM file for which you would like to calculate the mapping statistics?
read file

check=$(dirname $file)

if [ "$check" != '.' ]
then
    samFile=${file}
else
    samFile=${ngmDir}/${file}
fi

echo What is the name of the sample?
read sample

###############################
# Add -s value to output file #
###############################
naVal='N/A'

if [ "${sVal,,}" != "${naVal,,}" ]
then
    sval=$(printf '%.1f' $sVal)
else
    sval=$sVal
fi

###############################
# Count total number of reads #
###############################
numRawReads=$(cut -d ' ' -f2 ${projDir}/demultiplexed_reads_count.txt | awk '{SUM+=$1} END {print SUM/2}')
echo numRawReads=${numRawReads}

############################################
# Convert SAM to BAM and remove duplicates #
############################################
samFname=${samFile/.sam/}

samtools view -@ 18 -m 128G -1 $samFile -T $refSeq -b | samtools sort -n -@ 18 -m 128G - -o ${samFname}_namesorted.bam
samtools fixmate -m ${samFname}_namesorted.bam ${samFname}_namesorted_fixmate.bam
samtools sort ${samFname}_namesorted_fixmate.bam -o ${samFname}_positionsorted_fixmate.bam

bamFile=${samFname}_positionsorted_dedup.bam
samtools markdup -r ${samFname}_positionsorted_fixmate.bam $bamFile
samtools index -@ 18 -m 128G $bamFile ${bamFile}.bai

###########################################################################################################
# Count number of mapped reads aligning to chromosomes only, scaffolds only, and whole reference sequence #
###########################################################################################################
printf -v mappedChr "%d" $(samtools view -@ 18 -m 128G -F 4 -c $bamFile -L $bedChr)
echo mappedChr=${mappedChr}

printf -v mappedScaf "%d" $(samtools view -@ 18 -m 128G -F 4 -c $bamFile -L $bedScaf)
echo mappedScaf=${mappedScaf}

printf -v mappedTotal "%d" $(samtools view -@ 18 -m 128G -F 4 -c $bamFile)
echo mappedTotal=${mappedTotal}

################################################
# Count number of unmapped reads; SAM flag = 4 #
################################################
unmappedChr=$(samtools view -@ 18 -m 128G -f 4 -c $bamFile -L $bedChr)
echo unmappedChr=${unmappedChr}

unmappedScaf=$(samtools view -@ 18 -m 128G -f 4 -c $bamFile -L $bedScaf)
echo unmappedScaf=${unmappedScaf}

########################################
# Calculate % of reads that are mapped #
########################################
mapFracPerChr=$(echo $mappedChr $mappedTotal | awk '{printf "%.1f", ($1/$2)*100}')
echo mapFracPerChr=${mapFracPerChr}

mapFracPerScaf=$(echo $mappedScaf $mappedTotal | awk '{printf "%.1f", ($1/$2)*100}')
echo mapFracPerScaf=${mapFracPerScaf}

mapFracTotal=$(echo $mappedTotal $numRawReads | awk '{printf "%.10f", ($1/$2)*100}')
echo mapFracTotal=${mapFracTotal}

########################################################
# Count number of reads that uniquely map to reference #
########################################################
uniqReadsChr=$(samtools view -@ 18 -m 128G -F 4 $bamFile -L $bedChr | grep -e "NH\:i\:1\s" | wc -l)
echo uniqReadsChr=${uniqReadsChr}

uniqReadsScaf=$(samtools view -@ 18 -m 128G -F 4 $bamFile -L $bedScaf | grep -e "NH\:i\:1\s" | wc -l)
echo uniqReadsScaf=${uniqReadsScaf}

###################################################################
# Count number of reads that map to multiple regions of reference #
###################################################################
multiReadsChr=$(samtools view -@ 18 -m 128G -F 4 $bamFile -L $bedChr | grep -v "NH\:i\:[0-1]\s" | wc -l)
echo multiReadsChr=${multiReadsChr}

multiReadsScaf=$(samtools view -@ 18 -m 128G -F 4 $bamFile -L $bedScaf | grep -v "NH\:i\:[0-1]\s" | wc -l)
echo multiReadsScaf=${multiReadsScaf}

########################################################################
# Extract mapping quality scores (0 = non-unique, >10 probably unique) #
########################################################################
## Extract values of mapping quality scores and write to file
mapqChr=($(samtools view -@ 18 -m 128G -F 4 $bamFile -L $bedChr | grep -v "NH\:i\:0\s" | cut -f5)) # number of reads with non-zero MAPQ

mapqChrFile=${samFname}_mapq_scores_chr.txt
printf "%s\n" "${mapqChr[@]}" > $mapqChrFile

mapqScaf=($(samtools view -@ 18 -m 128G -F 4 $bamFile -L $bedScaf | grep -v "NH\:i\:0\s" | cut -f5))

mapqScafFile=${samFname}_mapq_scores_scaf.txt
printf "%s\n" "${mapqScaf[@]}" > $mapqScafFile

############################################
# Calculate average mapping quality scores #
############################################
mapqSumChr=$(paste -s -d+ $mapqChrFile | bc -l)
mapqNumChr=$(grep -cve "^\s*$" $mapqChrFile | bc -l)

printf -v mapqAvgChr "%.1f" $(( $mapqSumChr/$mapqNumChr ))
echo mapqAvgChr=${mapqAvgChr}

mapqSumScaf=$(paste -s -d+ $mapqScafFile | bc -l)
mapqNumScaf=$(grep -cve "^\s*$" $mapqScafFile | bc -l)

printf -v mapqAvgScaf "%.1f" $(( $mapqSumScaf/$mapqNumScaf ))
echo mapqAvgScaf=${mapqAvgScaf}

#####################################################
# Calculate breadth of coverage of reference genome #
#####################################################
## http://www.metagenomics.wiki/tools/samtools/breadth-of-coverage
## Get total number of bases sequenced at depth of 5X or higher
samtools depth ${samFname}_positionsorted_dedup.bam > ${samFname}_positionsorted_dedup_depth.txt

covFile=(${samFname}_positionsorted_dedup_depth.txt)
        
printf -v numBases5XcovChr "%d" $(cat $covFile | grep -E $chrNamesSearch | awk '$3>=5' | wc -l)
echo numBases5XcovChr=${numBases5XcovChr}

printf -v numBases5XcovScaf "%d" $(cat $covFile | grep -Ev $chrNamesSearch | awk '$3>=5' | wc -l)
echo numBases5XcovScaf=${numBases5XcovScaf}

printf -v numBases5XcovTotal "%d" $(cat $covFile | awk '$3>=5' | wc -l)
echo numBases5XcovTotal=${numBases5XcovTotal}

#########################################
# Calculate coverage breadth in percent #
#########################################
per5XcovChr=$(echo $numBases5XcovChr $refLenChr | awk -v v=$subsetPer '{printf "%.10f", ($1/($2*v))*100}')
echo per5XcovChr=${per5XcovChr} # % of mapped bases with >=5x coverage as proportion of total number of bases assigned to chromosomes in reference

per5XcovScaf=$(echo $numBases5XcovScaf $refLenScaf | awk -v v=$subsetPer '{printf "%.10f", ($1/($2*v))*100}')
echo per5XcovScaf=${per5XcovScaf}

per5XcovTotal=$(echo $numBases5XcovTotal $refLen | awk -v v=$subsetPer '{printf "%.10f", ($1/($2*v))*100}')
echo per5XcovTotal=${per5XcovTotal}

############################
# Add stats to output file #
############################
## Collect stats for each sample
row="${dataSet},${sample},${subsetPer},${sval},${minIdent},${numRawReads},${numReadsMerged},${mappedChr},${mappedScaf},${mappedTotal},${unmappedChr},${unmappedScaf},${mapFracPerChr},${mapFracPerScaf},${mapFracTotal},${uniqReadsChr},${uniqReadsScaf},${multiReadsChr},${multiReadsScaf},${mapqAvgChr},${mapqAvgScaf},${numBases5XcovChr},${numBases5XcovScaf},${numBases5XcovTotal},${per5XcovChr},${per5XcovScaf},${per5XcovTotal}"

############################
## Append stats to new row #
############################
echo ${row//,/$'\t'} >> $outTable

#################################################
# Write message after calculations are complete #
#################################################
echo -e "\033[1mcomputation of mapping stats complete: file=${samFile}\033[0m"