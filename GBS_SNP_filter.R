library(tidyverse)
parameters <- read.table("GBS_SNP_filter.txt")
basename <- gsub(".vcf","",parameters[1,1])
filelist <- list.files()

if (!((paste(basename,".biallelic.vcf",sep="")) %in% filelist)) {#1A  what to do if biallelic doesn't exist
  temp <- read_tsv("temp",col_names=TRUE)
  check <- temp %>% mutate(commacount=str_count(ALT,",")) %>% filter(., commacount==0) %>% 
  
  
  
 } #1B  
