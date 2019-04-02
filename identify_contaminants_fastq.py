#!/anaconda3/bin/python3.6

# Description: Use this script to examine FastQC reports and assess whether there are contaminants in your sequence data. Sequence files containing the forward and reverse reads must be indicated with "R1" and "R2" in the file names, respectively.
#
# Outputs: This script generates an HTML table listing potential sequence contaminants identified by FastQC.
#
# Software requirements: The following Python module is required - bs4 (BeautifulSoup).
#
# Usage: python identify_contaminants_fastq.py [full path name of FastQC directory] [suffix to append to output file name (e.g., unique identifier that specifies the project or analysis)]

##################
# Import modules #
##################
import codecs, glob, os, sys
from bs4 import BeautifulSoup

####################################################
# Pass FastQC outputs directory path to a variable #
####################################################
fastqcDir = sys.argv[1]

##############################################################
# Pass suffix to append to name of output file to a variable #
##############################################################
outSuffix = sys.argv[2]

#####################################
# Specify pipeline output directory #
#####################################
resultsDir = os.path.dirname(fastqcDir)

####################################
# Create list of FastQC html files #
####################################
reports = sorted(glob.glob("%s/*fastqc.html" % fastqcDir))

#################################################################
# Create sets to collect sequences of forward and reverse reads #
#################################################################
seqListR1 = set("")
seqListR2 = set("")

##############################################
# Find sequence contaminants and add to sets #
##############################################
for f in reports:
    page = codecs.open(f, "r", "utf-8")
    
    soup = BeautifulSoup(page, "html.parser")
    
    seqModule = soup.find(id = "M9")
    seqModuleHtml = seqModule.next_sibling
    seqModuleTable = seqModuleHtml.find_all("td")
    
    numSeqs = int(len(seqModuleTable)/4)
    indicesSeqs = [x*4 - 4 for x in range(1, numSeqs + 1)]
    indicesSeqNames = [x*4 - 1 for x in range(1, numSeqs + 1)]
    
    seqs = [str(seqModuleTable[i]).split(" (")[0].replace("<td>", "").replace("</td>", "").strip() for i in indicesSeqs]
    seqNames = [str(seqModuleTable[j]).split(" (")[0].replace("<td>", "").replace("</td>", "").strip() for j in indicesSeqNames]
    
    for k in range(0, len(seqs)):
        if "R1" in f:
            seqListR1.add("%s\t%s" % (seqNames[k], seqs[k]))
        
        elif "R2" in f:
            seqListR2.add("%s\t%s" % (seqNames[k], seqs[k]))

######################################################
# Write sequences of potential contaminants to files #
######################################################
if (seqListR1 != set() and seqListR2 != set()):
    adaptDict = {"R1": ("\n").join(sorted(list(seqListR1))), "R2": ("\n").join(sorted(list(seqListR2)))}
    
    for l in list(adaptDict.keys()):
        outtext = adaptDict[l]
        print("%s contaminants\n" % l, outtext, "\n")
        
        out = open("%s/fastqc_contaminants_%s_%s.txt" % (resultsDir, l, outSuffix), "w")
        out.write(outtext)
        out.close()

print("Sequences of potential contaminants identified by FastQC are written to files in directory '%s'" % resultsDir)
