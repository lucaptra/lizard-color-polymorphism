# Description: Use this script to generate an HTML table that summarizes assembly quality statistics calculated by Quast.
#
# Usage: Rscript summarize_quast_reports.R [species code] [full path name of Quast output directory]

##################
# Pass arguments #
##################
args <- commandArgs(trailingOnly = TRUE) # pass argument to script; args[1] = directory containing Quast outputs

###################
# Specify species #
###################
species <- args[1] # code used to identify species in sequencing data

#################################################################
# Get paths of Quast reports for each assembly for each species #
#################################################################
quastDir <- args[2] # full path of directory containing outputs from Quast

quastFiles <- Sys.glob(paste0(quastDir, "/summary_table_quast*", species, "*.txt"))

##################
# Load libraries #
##################
library(kableExtra)
library(rmarkdown)

###############################################
# Generate table summarizing Quast statistics #
###############################################
outTable <- data.frame() # create empty data frame

for (f in quastFiles) {
    t <- unname(read.table(f, header = FALSE, sep = " ", stringsAsFactors = FALSE, colClasses = c("character","integer", "integer", "character", "numeric", "numeric", "integer", "numeric")))
    names(t) <- c("assembly", "# contigs", "longest contig (bp)", "total length (bp)", "genome fraction (%)", "GC (%)", "N50", "duplication ratio")
    
    for (i in c(1:dim(t)[1])) {
        t[i, 5] <- formatC(t[i, 5], format = "e", digits = 2)
        
        outTable <- rbind(outTable, cbind(f, t[i, ]))
    }
}

names(outTable) <- c("file name", "assembly", "# contigs", "longest contig (bp)", "total length (bp)", "genome fraction (%)", "GC (%)", "N50", "duplication ratio") # name columns

htmlTable <- kable(outTable, caption = paste0("**Summary of _de novo_ assemblies for species ", species, "**"), format.args = list(big.mark = ","), align = "c") %>% kable_styling(bootstrap_options = "striped", full_width = F) # create kable table
save_kable(htmlTable, file = paste0(quastDir, "/quast_summary_table_", species, ".html"), self_contained = TRUE) # save table to HTML file

unlink(paste0(quastDir, "/lib"), recursive = TRUE) # delete folder created by save_kable function

print(paste0("HTML table summarizing assembly statistics for species ", species, " are written to directory '", quastDir, "'")) # print message after all tasks are complete