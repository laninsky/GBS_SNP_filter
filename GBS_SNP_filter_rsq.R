library(dplyr)
library(readr)
library(stringr)
#library(tidyverse) # I've had issues loading library(tidyverse) and R crashing using the sbatch syste
parameters <- read.table("GBS_SNP_filter.txt",header=FALSE,stringsAsFactors=FALSE)
basename <- gsub(".vcf","",parameters[1,1])
filelist <- list.files()
logfilename <- paste(basename,".log",sep="")

headerrows <- read_tsv("header_row.txt",col_names=FALSE)
numberofheaders <- dim(headerrows)[1]
temp <- read_tsv((paste(basename,".",parameters[2,1],"_",parameters[3,1],".vcf",sep="")),col_names=TRUE,skip=numberofheaders)
origcolnumber <- dim(temp)[2]
popmap <- read.table("popmap.txt",header=FALSE,stringsAsFactors=FALSE)
popnames <- unique(popmap[,2])

SNP_record <- NULL

for (k in 1:length(popnames)) {
  tempk <- read_table((paste(basename,".",parameters[2,1],"_",parameters[3,1],".",popnames[k],".pop.ld",sep="")),col_names=TRUE)
  tempk <- select(tempk,c(SNP_A,SNP_B,R2))
  tempk <- mutate_at(tempk,vars(SNP_A,SNP_B),funs(gsub(":.*","", . )))
  if(is.null(SNP_record)) {
    SNP_record <- tempk    
  }  else {
    SNP_record <- full_join(SNP_record,tempk,by=c("SNP_A","SNP_B"))
  } 
}
names(SNP_record) <- c("Locus_1","Locus_2",popnames)

as.numeric(parameters[6,1])
