# Description: Use this script to create an HTML table and plot a histogram of the number of fragments that are sequenced for a set of samples.
#
# Usage: Rscript calc_stats_read_counts.R [full path name of tab-delimited text file containing read counts]

##################
# Pass arguments #
##################
args <- commandArgs(trailingOnly = TRUE) # pass in argument from the command line; args[1] = full path to the tab-delimited file

#################################
# Set pipeline output directory #
#################################
resultsDir <- dirname(args[1]) # get output directory

###################################
# Import text file of read counts #
###################################
t <- read.table(args[1], header = TRUE, sep = "\t", colClasses = c("character", "numeric"), stringsAsFactors = FALSE) # import count data

library(textclean)
t[ , 1] <- mgsub(t[ , 1], c("_R.", ".fastq.*?"), c("", ""), fixed = FALSE) # extract sample names and substitute values in column 1

names(t) <- c("sample_name", "number_sequenced_fragments") # set column names

t <- aggregate(number_sequenced_fragments ~ sample_name, data = t, FUN = function(x) sum(x)/2) # calculate number of fragments [(# F reads + # R reads)/2]

##########################################
# Generate table summarizing read counts #
##########################################
outFname <- strsplit(basename(args[1]), ".", fixed = TRUE)[[1]][1] # extract name of output file

library(kableExtra) # load library to export pretty R table

htmlTable <- kable(t, format = "html", digits = 0, row.names = FALSE, col.names = c("sample_name", "#_sequenced_fragments"), align = "c", caption = "No. of sequenced fragments per sample", format.args = list(big.mark = ",")) %>% kable_styling(bootstrap_options = "striped", full_width = F) # organize data into table

save_kable(htmlTable, file = paste0(resultsDir, "/", outFname, "_table.html"), self_contained = TRUE) # save table to HTML file

unlist(paste0(resultsDir, "/lib"), recursive = TRUE) # delete folder created by save_kable function

############################################
# Plot histogram of per-sample read counts #
############################################
png(file = paste0(resultsDir, "/", outFname, "_hist_plot.png"), width = 11, height = 8.5, units = "in", pointsize = 24, res = 300) # write histogram to PNG file
par(mar = c(2, 4, 1, 1))
hist(t$number_sequenced_fragments, breaks = 15, main = "", xlab = expression(bold("no. of sequenced fragments per sample")), ylab = expression(bold("frequency")))
dev.off()

print(paste0("Table and histogram of the number of sequenced fragments per sample saved to directory '", resultsDir, "'!")) # print message after all tasks are complete