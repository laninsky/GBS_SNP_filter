#GBS_SNP_filter v1.16
# Loading in required packages
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('readr')) install.packages('readr'); library('readr')
if (!require('stringr')) install.packages('stringr'); library('stringr')
if (!require('tibble')) install.packages('tibble'); library('tibble')
#library(tidyverse) # I've had issues loading library(tidyverse) and R crashing using the sbatch system

mapfile <- read_tsv(list.files(pattern="*.map$"),col_names=FALSE)

if (is.numeric(mapfile$X1)) {
  mapfile$X1 <- paste("GBS_SNP_filter_",mapfile$X1,sep="")
}

write_tsv(mapfile,list.files(pattern="*.map$"),col_names=FALSE)

print("GBS_SNP_filter_chrom_modifier.R has the following warnings():")
warnings()
print("If no warnings printed, none thrown by GBS_SNP_filter_chrom_modifier.R")
