library(tidyverse)
parameters <- read.table("GBS_SNP_filter.txt")
basename <- gsub(".vcf","",parameters[1,1])
filelist <- list.files()
logfilename <- paste(basename,".log",sep="")

if (!((paste(basename,".biallelic.vcf",sep="")) %in% filelist)) {#1A  what to do if biallelic doesn't exist
  temp <- read_tsv("temp",col_names=TRUE)
  write(format(Sys.time(),usetz = TRUE),logfilename,append=TRUE)
  write(paste(basename,".vcf has the following number of SNPs:",sep=""),logfilename,append=TRUE)
  write((dim(temp)[1]-1),logfilename,append=TRUE)  
  temp <- temp %>% mutate(commacount=str_count(ALT,",")) %>% filter(., commacount==0) %>% select(.,-commacount)
  write(format(Sys.time(),usetz = TRUE),logfilename,append=TRUE)
  write(paste(basename,".biallelic.vcf has the following number of SNPs:",sep=""),logfilename,append=TRUE)
  write((dim(temp)[1]-1),logfilename,append=TRUE)  
  write_delim(temp,(paste(basename,".biallelic.vcf",sep="")),delim="\t",append=TRUE,col_names=TRUE)  
 } #1B  
