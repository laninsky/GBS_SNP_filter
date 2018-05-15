library(tidyverse)
parameters <- read.table("GBS_SNP_filter.txt")
basename <- gsub(".vcf","",parameters[1,1])
filelist <- list.files()
logfilename <- paste(basename,".log",sep="")

if (!((paste(basename,".biallelic.vcf",sep="")) %in% filelist)) {#1A  what to do if biallelic doesn't exist
  temp <- read_tsv("temp",col_names=TRUE)
  write(format(Sys.time(),usetz = TRUE),logfilename,append=TRUE)
  write(paste(basename,".vcf has the following number of SNPs:",sep=""),logfilename,append=TRUE)
  write((dim(temp)[1]),logfilename,append=TRUE)  
  temp <- temp %>% mutate(commacount=str_count(ALT,",")) %>% filter(., commacount==0) %>% select(.,-commacount)
  write(format(Sys.time(),usetz = TRUE),logfilename,append=TRUE)
  write(paste(basename,".biallelic.vcf has the following number of SNPs:",sep=""),logfilename,append=TRUE)
  write((dim(temp)[1]),logfilename,append=TRUE) 
  headerrows <- read_tsv("header_row.txt",col_names=FALSE)
  write.table(headerrows,(paste(basename,".biallelic.vcf",sep="")),quote=FALSE,row.names=FALSE,col.names=FALSE)  
  write_delim(temp,(paste(basename,".biallelic.vcf",sep="")),delim="\t",append=TRUE,col_names=TRUE)  
} else { #1AB
  if (!((paste(basename,".oneSNP.vcf",sep="")) %in% filelist)) {#2A reading in *biallelic.vcf if *oneSNP.vcf doesn't exist
    headerrows <- read_tsv("header_row.txt",col_names=FALSE)
    numberofheaders <- dim(headerrows)[1]
    temp <- read_tsv((paste(basename,".biallelic.vcf",sep="")),col_names=TRUE,skip=numberofheaders)
  } #2B  
} #1B  

if (!((paste(basename,".oneSNP.vcf",sep="")) %in% filelist)) {#3A: if oneSNP.vcf doesn't exist, filtering for SNPs with the greatest coverage across individuals
  duplicatedloci <- unique(temp$`#CHROM`[which(duplicated(temp$`#CHROM`)==TRUE)])
  write(format(Sys.time(),usetz = TRUE),logfilename,append=TRUE)
  write("The following number of loci have more than one SNP:",logfilename,append=TRUE) 
  write(length(duplicatedloci),logfilename,append=TRUE)
  notduplicated <- temp %>% filter(., (!(`#CHROM` %in% duplicatedloci)))
  notduplicated <- notduplicated %>%  mutate_at(vars(10:dim(temp)[2]), .funs = funs(cov = gsub(":.*","",gsub("^.*?:","", . )))) 
  notduplicated <- notduplicated %>%  mutate_at(vars((dim(temp)[2]+1):(dim(notduplicated)[2])),funs(as.numeric)) 
  duplicated <- temp %>% filter(., (`#CHROM` %in% duplicatedloci))
  duplicated <- duplicated %>% mutate_at(vars(10:dim(temp)[2]), .funs = funs(cov = gsub(":.*","",gsub("^.*?:","", . )))) 
  duplicated <- duplicated %>%  mutate_at(vars((dim(temp)[2]+1):(dim(duplicated)[2])),funs(as.numeric))
  duplicated_reduced <- NULL
  for (i in duplicatedloci) {
    temptemp <- duplicated %>% filter(., (`#CHROM` %in% i))
    zerocounts <- unlist(lapply(1:(dim(temptemp)[1]),function(x){sum(temptemp[x,(dim(temp)[2]+1):(dim(temptemp)[2])]==0)})) 
    temptemp <- temptemp[(which(zerocounts==min(zerocounts))),]    
    covcounts <- rowSums(temptemp[,(dim(temp)[2]+1):(dim(temptemp)[2])])
    temptemp <- temptemp[(which(covcounts==max(covcounts))),]
    if (dim(temptemp)[1]>1) { #4A: if we haven't narrowed down to one SNP for this locus
      randomSNP <- sample(1:dim(temptemp)[1], 1, replace=FALSE)
      duplicated_reduced <- bind_rows(duplicated_reduced, temptemp[randomSNP,])
    } else { #4AB: if we HAVE narrowed to one SNP for this locus
      duplicated_reduced <- bind_rows(duplicated_reduced, temptemp)
    } #4B 
    print(paste("Up to ",dim(duplicated_reduced)[1]," out of ",length(duplicatedloci)," loci with multiple SNPs",sep=""))
  }
  origcolnumber <- dim(temp)[2]
  temp <- bind_rows(notduplicated,duplicated_reduced)
  write(format(Sys.time(),usetz = TRUE),logfilename,append=TRUE)
  write(paste(basename,".oneSNP.vcf has the following number of SNPs:",sep=""),logfilename,append=TRUE)
  write((dim(temp)[1]),logfilename,append=TRUE) 
  write.table(headerrows,(paste(basename,".oneSNP.vcf",sep="")),quote=FALSE,row.names=FALSE,col.names=FALSE)  
  write_delim(temp[,1:origcolnumber],(paste(basename,".oneSNP.vcf",sep="")),delim="\t",append=TRUE,col_names=TRUE)  
} else { #3AB if oneSNP.vcf does exist
    headerrows <- read_tsv("header_row.txt",col_names=FALSE)
    numberofheaders <- dim(headerrows)[1]
    temp <- read_tsv((paste(basename,".oneSNP.vcf",sep="")),col_names=TRUE,skip=numberofheaders)
    origcolnumber <- dim(temp)[2]
    temp <- temp %>% mutate_at(vars(10:dim(temp)[2]), .funs = funs(cov = gsub(":.*","",gsub("^.*?:","", . )))) 
    temp <- temp %>%  mutate_at(vars((origcolnumber+1):(dim(temp)[2])),funs(as.numeric))
} #3B

if (!((paste(basename,".rsq",sep="")) %in% filelist)) { #4A: if *.rsq doesn't exist, creating this
  popmap <- read_tsv("popmap.txt",col_names=FALSE)
  
  # need to read in popmap and give rsq by population between SNPs
  # e.g. locus_1 locus_2 pop1 pop2 pop3 etc
} #4B  
  
# after all of this preamble, the actual code for doing   
  
