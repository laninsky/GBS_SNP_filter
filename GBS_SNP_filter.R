library(tidyverse)
parameters <- read.table("GBS_SNP_filter.txt")
basename <- gsub(".vcf","",parameters[1,1])
filelist <- list.files()

if (!((paste(basename,".biallelic.vcf",sep="")) %in% filelist)) {#1A  what to do if biallelic doesn't exist
  temp <- read_table("temp",col_names=TRUE)
  temp <- separate(temp,col=1,sep="\t",remove=TRUE)
  temp <- as_tibble(read.table("temp",header=TRUE,stringsAsFactors=FALSE))
} #1B  
