#!/bin/bash

# Description: Use this script to run Quast to evaluate the quality of de novo genome assemblies.
#
# Outputs: This script generates text and HTML reports summarizing various assembly statistics calculated by Quast.
#
# Software requirements: Quast
#
# Usage: bash assess_denovo_assembly_quality_quast.sh [species code] [name of de novo assembler] [full path of project directory] [name of reference sequence]

###################
# Specify species #
###################
species=JAR # code used to identify species in sequencing data

#############################
# Specify de novo assembler #
#############################
assembler=mccortex

#############################
# Specify input directories #
#############################
projDir=/scratch/ltran20/lizards # directory containing all project directories

mccpipdir=${projDir}/mccortex_uncorr_${species}/ # directory containing outputs from McCortex

#######################################
# Specify directory for Quast outputs #
#######################################
quastDir=${projDir}/quast_results
mkdir -p $quastDir # create directory if it doesn't already exist

############################
# Specify reference genome #
############################
reference=(/home/ltran20/lizards/assemblies/*.fa) # full path of reference FASTA

###########################################
# Delete any previous Quast summary table #
###########################################
outTable=${quastDir}/summary_table_quast_${assembler,,}_${species}_full.txt

if [ ! -f "$outTable" ]
then
    :
else
    rm $outTable
fi

##############################
# Run Quast on each assembly #
##############################
if [ "${assembler,,}" = "mccortex" ]
then
    # Get McCortex contig assemblies
    declare -a contigFiles=(${mccpipdir}*rmdup*.fa.gz)
    
    # For each assembly
    for fasta in "${contigFiles[@]}"
    do        
        fname=$(basename $fasta)
        sample=${fname/.pop.rmdup.contigs.fa.gz/}
        
        rm -Rf ${quastDir}/quast_${assembler,,}_${sample}_full
        
        python quast.py $fasta -o ${quastDir}/quast_${assembler,,}_${sample}_full -R $reference --min-contig 50 -t 18 -l ${sample}_full --eukaryote --plots-format png --space-efficient # run Quast
        
        # Extract assembly stats from Quast report
        quastReport=(${quastDir}/quast_${assembler,,}_${sample}_full_report.tsv)
        
        row=()
        
        for pattern in 'Assembly' '# contigs (>= 0 bp)' 'Largest contig' 'Total length (>= 0 bp)' 'Genome fraction (%)' '^GC (%)' 'N50' 'Duplication ratio'
        do
            row+=("$(grep -e "${pattern}" $quastReport | cut -d$'\t' -f 2)")
        done
        
        echo ${row[@]} >> $outTable
    done
fi