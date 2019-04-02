# Description: Use this script to create an HTML table that summarizes quality control statistics calculated by FastQC/MultiQC.
#
# Usage: Rscript summarize_fastqc_reports.R [full path name of FastQC output directory]

##################
# Pass arguments #
##################
args <- commandArgs(trailingOnly = TRUE) # pass argument to script; args[1] = directory containing FastQC reports and MultiQC reports for both forward and reverse reads

################################
# Get paths of MultiQC reports #
################################
resultsDir <- dirname(args[1])

multiqcFiles <- Sys.glob(paste0(args[1], "/*/multiqc_fastqc.txt")) # get path names of MultiQC reports for reads of either orientation

##################
# Load libraries #
##################
library(kableExtra)
library(lettercase)
library(plyr)
library(textclean)

##############################################
# Generate table summarizing MultiQC reports #
##############################################
for (f in multiqcFiles) {
    outFname <- basename(dirname(f)) # extract name of output file
    
    readOrient <- gsub("multiqc_data_", "", basename(dirname(f))) # get orientation of reads
    
    t <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE) # import MultiQC report
    
    names(t) <- gsub("\\.", "_", str_lowercase(names(t))) # edit column names to comply with R conventions
    names(t)[8] <- "%_GC"
    
    t[ , 1] <- mgsub(t[ , 1], c("_R.", ".fastq.*?"), c("", ""), fixed = FALSE) # extract sample names and substitute values in column 1
    
    formatCmd <- c() # set font color to red for entries in which quality control test failed or is close to failing
    for (n in names(t)[11:21]) {
        formatCmd <- c(formatCmd, paste0(n, " = cell_spec(", n, ", 'html', color = ifelse(", n, " == 'warn', 'red', 'black'))"))
    }
    
    kableCmd <- paste0("t %>% mutate(", paste(formatCmd, collapse = ", "), ") %>% kable(format = 'html', digits = 2, row.names = FALSE, align = 'c', caption = 'Summary of quality control checks on fastq data (", readOrient, " reads)', format.args = list(big.mark = ','), escape = F) %>% kable_styling(bootstrap_options = 'striped', full_width = F)") # command to format HTML table
    
    htmlTable <- eval(parse(text = kableCmd))
    
    save_kable(htmlTable, file = paste0(resultsDir, "/", outFname, "_", readOrient, "_table.html"), self_contained = TRUE) # save table to HTML file
    
    unlist(paste0(resultsDir, "/lib"), recursive = TRUE) # delete folder created by save_kable function
}

print(paste0("HTML tables summarizing quality control statistics are written to directory '", resultsDir, "'")) # print message after all tasks are complete