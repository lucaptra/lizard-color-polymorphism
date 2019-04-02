# Description: Use this script to generate various plots and statistics that can be used to assess the quality of the adapter removal and/or read trimming step(s).
#
# Usage: Rscript calc_stats_plots_cleaned_reads.R [full path name of project directory] [name of file containing number of mapped reads per sample] [name of directory containing files of read lengths]

##################
# Pass arguments #
##################
args <- commandArgs(trailingOnly = TRUE) # arguments to pass to script

#############################
# Specify directories #
#############################
projDir <- args[1]

resultsDir <- paste0(projDir, "/processing_results")

####################################################
# Import data on number of mapped reads per sample #
####################################################
t1 <- read.table(paste0(resultsDir, "/", args[2]), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
names(t1) <- c("sample", "no._mapped_reads")

########################################################
# Get list of files containing per-sample read lengths #
########################################################
files <- list.files(path = paste0(resultsDir, "/", args[3]), pattern = "*_read_lengths.txt", full.names = TRUE)

#################################################
# Calculate average read length for each sample #
#################################################
avgReadLenDF <- data.frame()

for (f in files) {
    sample = gsub("_read_lengths.txt", "", basename(f))
    
    t2 <- read.table(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE, colClasses = "integer")
    
    lengths <- as.numeric(unlist(t2))
    
    avgReadLen <- sum(lengths)/length(lengths)
    avgReadLenDF <- rbind(avgReadLenDF, data.frame(sample = sample, avg_read_length_bp = avgReadLen))
}

write.table(avgReadLenDF, paste0(resultsDir, "/average_read_lengths_samples.txt"), quote = FALSE, sep = "\t", row.names = FALSE)

################################################################################################
# Merge data frames containing per-sample lengths of cleaned reads and numbers of mapped reads #
################################################################################################
finalDF <- Reduce(function(x, y) merge(x, y, by = "sample", all.y = TRUE), list(t1, avgReadLenDF))

######################################################################
# Plot scatterplot of cleaned read lengths v. number of mapped reads #
######################################################################
png(paste0(resultsDir, "/cleaned_read_lengths_mapping_scatterplots_samples.png"), width = 11, height = 8.5, units = "in", res = 300)

plot(finalDF$no._mapped_reads ~ finalDF$avg_read_length_bp, pch = 20, main = "", xlab = expression(bold("mean cleaned read length per sample")), ylab = expression(bold("no. of mapped reads per sample")))
with(finalDF, text(no._mapped_reads ~ avg_read_length_bp, labels = finalDF$sample, pos = 2))

dev.off()

#####################################################################
# Regress per-sample number of mapped reads on cleaned read lengths #
#####################################################################
lmSummary <- summary(lm(finalDF$no._mapped_reads ~ finalDF$avg_read_length_bp))

print(paste0("slope = ", lmSummary$coefficients[2, 1]))
print(paste0("*p*-value = ", signif(lmSummary$coefficients[2, 4], digits = 2)))