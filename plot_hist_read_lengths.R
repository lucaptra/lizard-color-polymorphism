# Description: Use this script to plot a histogram of the lengths of reads in fastq/fastq.gz files in a specified directory.
#
# Usage: Rscript plot_hist_read_lengths.R [full path name of directory containing per-sample read length files]

##################
# Pass arguments #
##################
args <- commandArgs(trailingOnly = TRUE) # path name of directory containing files summarizing read lengths

#####################################
# Specify pipeline output directory #
#####################################
resultsDir <- args[1]

#############################################
# Get list of files containing read lengths #
#############################################
files <- list.files(path = resultsDir, pattern = "*_read_lengths.txt", full.names = TRUE)

##################################
# Plot histogram of read lengths #
##################################
png(paste0(resultsDir, "/hist_trimmed_read_lengths_per_sample.png"), width = 11, height = 8.5, units = "in", res = 300)
par(mfrow = c(3, 6), oma = c(1, 1.3, 0, 0))

for (f in files) {
    t <- read.table(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE, colClasses = "integer")
    hist(as.numeric(unlist(t)), xlim = c(0, 200), ylim = c(0, 70000000), main = i, xlab = "", ylab = "")
}

mtext("read lengths", side = 1, outer = TRUE, font = 2)
mtext("frequency", side = 2, outer = TRUE, font = 2)

dev.off()