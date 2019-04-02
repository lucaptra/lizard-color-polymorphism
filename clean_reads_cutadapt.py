#!/anaconda3/bin/python3.6

# Description: Use this script to remove contaminant adapter sequences using CutAdapt from your sequence data.
#
# Outputs: This script generates fastq files containing cleaned reads. The orientation of your reads must be indicated with "R1" or "R2" in the file names.
#
# Software requirements: CutAdapt (https://cutadapt.readthedocs.io/en/stable/guide.html#removing-adapters) must be installed. The following Python modules are also required - Biopython, wget.
#
# Usage: python clean_reads_cutadapt.py [full path name of directory containing raw/demultiplexed reads] [name of directory for CutAdapt outputs] [prefix of the adapters you want removed]

##################
# Import modules #
##################
import glob, heapq, os, shutil, subprocess as sub, sys, wget
from Bio import Seq, SeqIO

#########################################################################################
# Pass path name of directory in which raw/demultiplexed reads are stored to a variable #
#########################################################################################
rawReadsDir = sys.argv[1]

#####################################
# Specify pipeline output directory #
#####################################
resultsDir = "%s/processing_results" % os.path.dirname(rawReadsDir)

##################################################################################
# Pass path name of directory in which cleaned reads will be saved to a variable #
##################################################################################
cleanDir = "%s/%s" % (resultsDir, sys.argv[2])

###########################################################################
# Provide string containing identifier to isolate adapter(s) to be tested #
###########################################################################
adaptPrefix = sys.argv[3]

#######################################################
# Create output directory if it doesn't already exist #
#######################################################
if os.path.exists(cleanDir) != True:
    os.mkdir(cleanDir)
else:
    shutil.rmtree(cleanDir)
    os.mkdir(cleanDir)

##########################################################
# Download FASTA file listing adapter sequences to check #
##########################################################
filename = "%s/sequencing_adapters.fa" % resultsDir

url = "https://raw.githubusercontent.com/NCBI-Hackathons/OnlineAdapterDatabase/master/datasources/sequencing_adapters.fa"
wget.download(url, filename)

#########################
# Get names of adapters #
#########################
getAdaptSeqNames = sub.Popen('grep ">%s" %s | sed -e "s/>//g"' % (adaptPrefix, filename), shell = True, stdout = sub.PIPE)
seqNames = getAdaptSeqNames.stdout.read().decode("utf-8").strip().split("\n")

#############################
# Extract adapter sequences #
#############################
seqsF = [] # create a list to store the forward sequences

for record in SeqIO.parse(filename, "fasta"): # create an iterable sequence object
    if record.id in seqNames:
        seqsF.append(str(record.seq)) # add the sequence to the list
    if adaptPrefix.casefold() == "Truseq".casefold():
        if record.id == "TruSeqUniversalAdapter":
            rc = str(Seq.reverse_complement(record.seq)) # find the reverse complement of the Truseq universal adapter

#############################################################
# Write string joining adapter sequences to run in CutAdapt #
#############################################################
seqsF.insert(0, "")

adaptStr = (" -a ").join(seqsF)
adaptStr + " -A %s" % rc

##########################################################################
# Determine substring to remove from file names to retrieve sample names #
##########################################################################
files = glob.glob("%s/*_R1_*fastq*" % rawReadsDir) # list of paths to the sequence files that will be processed with CutAdapt
filenames = sorted([os.path.basename(x) for x in files]) # extract base names from the file paths

nameLens = [len(i) for i in filenames] # list of lengths of the file base names

substrList = [] # empty list to store substrings

for j in filenames: # determine all possible substrings of each base name up to the maximum length, which is set to the minimum value of the base name lengths
    for k in range(0, min(nameLens)):
        substrList.append(j[-k:])

substrDict = {} # create dictionary storing the number of occurrences of each possible substring

for l in substrList:
    substrDict[l] = substrList.count(l) # key = substring, value = number of occurrences

maxSubstrOccur = heapq.nlargest(1, list(set(list(substrDict.values())))) # calculate the largest and second largest number of occurrences of a substring

substrReplaceDict = {} # dictionary to store final set of substrings and their replacements

for m in maxSubstrOccur: # get substrings that occurred the most number of times
    maxSubstrDict = { key:value for key, value in substrDict.items() if value == m } # keep only those substrings in the dictionary that occur the most number of times
    maxSubstrKeys = list(maxSubstrDict.keys()) # get the substrings of the subsetted dictionary
    
    for n in range(0, len(maxSubstrKeys)):
        substrReplaceDict[maxSubstrKeys[n]] = "" # append to final dictionary of substrings

############################
# Print substrings to file #
############################
substrR1 = "\n".join(list(substrReplaceDict))

outFile = open("%s/filename_substrings_to_remove.txt" % resultsDir, "w")
outFile.write("\n".join(substrR1))
outFile.close()

###################################################
# Clean reads of adapter sequences using CutAdapt #
###################################################
for f1 in files:
    for o, p in substrReplaceDict.items():
        sample = os.path.basename(f1).replace(o, p)
    
    f2 = glob.glob("%s/%s*_R2_*fastq*" % (rawReadsDir, sample))[0]
    
    cmd = "cutadapt%s --times 26 --minimum-length 20 -q 10 -o %s/%s_R1_cutadapt.fastq.gz -p %s/%s_R2_cutadapt.fastq.gz %s %s --cores=18" % (adaptStr, cleanDir, sample, cleanDir, sample, f1, f2)
    p = sub.Popen(cmd, shell = True)
    p.communicate()
    
    print("Adapter sequences removed from sequence files of sample %s." % sample)
