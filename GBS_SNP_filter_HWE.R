library(dplyr)
library(readr)
library(stringr)
#library(tidyverse) # I've had issues loading library(tidyverse) and R crashing using the sbatch syste
parameters <- read.table("GBS_SNP_filter.txt",header=FALSE,stringsAsFactors=FALSE)
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
    temp <- temp %>% mutate_at(vars((origcolnumber+1):(dim(temp)[2])),funs(as.numeric))
} #3B

if (!((paste(basename,".",parameters[2,1],"_",parameters[3,1],".vcf",sep="")) %in% filelist)) { #4A: if dataset hasn't previously been filtered with the same parameters
  zerocounts <- unlist(lapply(1:(dim(temp)[1]),function(x){sum((temp[x,(origcolnumber+1):(dim(temp)[2])])==0)}))
  keeprows <- which(zerocounts<=((1-as.numeric(parameters[2,1]))*((dim(temp)[2])-origcolnumber)))
  temp <- temp[keeprows,]
  write(format(Sys.time(),usetz = TRUE),logfilename,append=TRUE)
  write(paste("Following filtering of SNPs not present in at least ",(as.numeric(parameters[2,1])*100),"% of samples, the following number of SNPs remained:",sep=""),logfilename,append=TRUE)
  write((dim(temp)[1]),logfilename,append=TRUE)   
  zerocounts <- unlist(lapply((origcolnumber+1):(dim(temp)[2]),function(x){sum((temp[1:(dim(temp)[1]),x])==0)}))
  removecols <- which(zerocounts>((as.numeric(parameters[3,1]))*(dim(temp)[1])))                  
  write(format(Sys.time(),usetz = TRUE),logfilename,append=TRUE)
  write(paste("The following samples have been removed because they have >",(as.numeric(parameters[3,1])*100), "% missing data:",sep=""),logfilename,append=TRUE)
  write((names(temp)[(removecols+9)]),logfilename,append=TRUE)
  temp <- select(temp,-one_of(names(temp)[c((removecols+9),(removecols+origcolnumber))]))
  origcolnumber <- origcolnumber-length(removecols) 
  write.table(headerrows,(paste(basename,".",parameters[2,1],"_",parameters[3,1],".vcf",sep="")),quote=FALSE,row.names=FALSE,col.names=FALSE)  
  write_delim(temp[,1:origcolnumber],(paste(basename,".",parameters[2,1],"_",parameters[3,1],".vcf",sep="")),delim="\t",append=TRUE,col_names=TRUE)    
  write(format(Sys.time(),usetz = TRUE),logfilename,append=TRUE)  
  write(paste("Following this filtering ",basename,".",parameters[2,1],"_",parameters[3,1],".vcf has been written out, containing ",(dim(temp)[1])," SNPs and ", (origcolnumber-9), " samples",sep=""),logfilename,append=TRUE)
} else { 
  headerrows <- read_tsv("header_row.txt",col_names=FALSE)
  numberofheaders <- dim(headerrows)[1]
  temp <- read_tsv((paste(basename,".",parameters[2,1],"_",parameters[3,1],".vcf",sep="")),col_names=TRUE,skip=numberofheaders)
  origcolnumber <- dim(temp)[2]
  temp <- temp %>% mutate_at(vars(10:dim(temp)[2]), .funs = funs(cov = gsub(":.*","",gsub("^.*?:","", . )))) 
  temp <- temp %>% mutate_at(vars((origcolnumber+1):(dim(temp)[2])),funs(as.numeric))
}  #4B

if (!((paste(basename,".",parameters[2,1],"_",parameters[3,1],".",parameters[4,1],"_",parameters[6,1],".HWE.vcf",sep="")) %in% filelist)) { #5A: If we haven't carried out HWE filtering for files with this combo of parameters yet
  if (!((paste(basename,".",parameters[2,1],"_",parameters[3,1],".",parameters[4,1],"_",parameters[6,1],".HWE",sep="")) %in% filelist)) { # 6A: If locus specific HWE values have not already been printed out  
    popmap <- read.table("popmap.txt",header=FALSE,stringsAsFactors=FALSE)
    popnames <- unique(popmap[,2])
    temptemp <- temp[,1:origcolnumber] 
    hwetable <- NULL
    for (k in 1:length(popnames)) { #8A: for each population
      print(paste("Up to ",k," out of ",length(popnames), " populations",sep="")) 
      temptemppop <- select(temptemp, which(names(temptemp) %in% (popmap[(which(popmap[,2]==popnames[k])),1])))
      temptemppop <- mutate_at(temptemppop,vars(1:dim(temptemppop)[2]),funs(gsub(":.*","", . )))
      hwetablepvalues <- unlist(lapply(1:(dim(temp)[1]),function(x){
        tempmatrix <- matrix(0,ncol=2,nrow=3)
        tempmatrix[1,1] <- length(which(temptemppop[x,]=="0/0"))
        tempmatrix[2,1] <- length(which((temptemppop[x,]=="0/1" | temptemppop[x,]=="1/0")))
        tempmatrix[3,1] <- length(which(temptemppop[x,]=="1/1"))
        tempmatrix[1,2] <- ((((2*tempmatrix[1,1])+tempmatrix[2,1])/(2*sum(tempmatrix[,1])))^2)*sum(tempmatrix[,1])
        tempmatrix[3,2] <- ((((2*tempmatrix[3,1])+tempmatrix[2,1])/(2*sum(tempmatrix[,1])))^2)*sum(tempmatrix[,1])
        tempmatrix[2,2] <- 2*(((2*tempmatrix[1,1])+tempmatrix[2,1])/(2*sum(tempmatrix[,1])))*(((2*tempmatrix[3,1])+tempmatrix[2,1])/(2*sum(tempmatrix[,1])))*sum(tempmatrix[,1])
        if (sum(tempmatrix[,1])==0) {
          temphwep <- "NaN"
        } else {  
          temphwep <- suppressWarnings(fisher.test(tempmatrix)$p.value)
        }       
        return(temphwep)
      }))
      hwetable <- cbind(hwetable,hwetablepvalues)
    } #8B  
    hwetable <- cbind(as.matrix(temp[,1]),hwetable)
    write.table(matrix(c("snp",popnames),nrow=1),(paste(basename,".",parameters[2,1],"_",parameters[3,1],".",parameters[4,1],"_",parameters[6,1],".HWE",sep="")),quote=FALSE,row.names=FALSE,col.names=FALSE) 
    write.table(hwetable,(paste(basename,".",parameters[2,1],"_",parameters[3,1],".",parameters[4,1],"_",parameters[6,1],".HWE",sep="")),quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
  } else { #6AB reading in existing HWE file  
    hwetable <- read.table(paste(basename,".",parameters[2,1],"_",parameters[3,1],".",parameters[4,1],"_",parameters[6,1],".HWE",sep=""),header=TRUE,stringsAsFactors=FALSE)
  }  
  
  outofhwe <- unlist(lapply(1:(dim(temp)[1]),function(x){
    sum(as.numeric(hwetable[x,2:(dim(hwetable)[2])])<as.numeric(parameters[4,1]))>=as.numeric(parameters[6,1]) 
  }))  
  
  hwetablebin <- hwetable[(which(outofhwe==TRUE)),]
  hwetablebin <- as.data.frame(hwetablebin)
  names(hwetablebin) <- c("snp",popnames)
  
  write(format(Sys.time(),usetz = TRUE),logfilename,append=TRUE)
  write(paste("The following loci (p-values given) will be removed as more than ",parameters[6,1]," populations had a HWE p-value of <",parameters[4,1],sep=""),logfilename,append=TRUE)
  write.table(hwetablebin,logfilename,append=TRUE,row.names=FALSE,col.names=TRUE,quote=FALSE)  
  temp <- temp %>% filter(., (!(`#CHROM` %in% (hwetablebin[,1]))))
  write.table(headerrows,(paste(basename,".",parameters[2,1],"_",parameters[3,1],".",parameters[4,1],"_",parameters[6,1],".HWE.vcf",sep="")),quote=FALSE,row.names=FALSE,col.names=FALSE)
  write_delim(temp[,1:origcolnumber],(paste(basename,".",parameters[2,1],"_",parameters[3,1],".",parameters[4,1],"_",parameters[6,1],".HWE.vcf",sep="")),delim="\t",append=TRUE,col_names=TRUE)    
  write(format(Sys.time(),usetz = TRUE),logfilename,append=TRUE)  
  write(paste("Following this filtering ",basename,".",parameters[2,1],"_",parameters[3,1],".",parameters[4,1],"_",parameters[6,1],".HWE.vcf has been written out, containing ",(dim(temp)[1])," SNPs and ", (origcolnumber-9), " samples",sep=""),logfilename,append=TRUE)   
} else {  #5AB
  headerrows <- read_tsv("header_row.txt",col_names=FALSE)
  numberofheaders <- dim(headerrows)[1]
  temp <- read_tsv((paste(basename,".",parameters[2,1],"_",parameters[3,1],".",parameters[4,1],"_",parameters[6,1],".HWE.vcf",sep="")),col_names=TRUE,skip=numberofheaders)
  origcolnumber <- dim(temp)[2]
  temp <- temp %>% mutate_at(vars(10:dim(temp)[2]), .funs = funs(cov = gsub(":.*","",gsub("^.*?:","", . )))) 
  temp <- temp %>% mutate_at(vars((origcolnumber+1):(dim(temp)[2])),funs(as.numeric))
} #5B


popmap <- read.table("popmap.txt",header=FALSE,stringsAsFactors=FALSE)
popnames <- unique(popmap[,2])
for (k in 1:length(popnames)) {
   tempK <- select(temp, c((1:10),which(names(temp) %in% (popmap[(which(popmap[,2]==popnames[k])),1]))))
   origcolnumber <- dim(tempK)[2]
   tempK <- mutate_at(tempK,vars(10:origcolnumber),.funs = funs(genotype = gsub(":.*","", . )))     
   tempK <- mutate(tempK, hom1 = rowSums(tempK[,(origcolnumber+1):(dim(tempK)[2])] == "1/1"))
   tempK <- mutate(tempK, het = ((rowSums(tempK[,(origcolnumber+1):(dim(tempK)[2]-1)] == "0/1")+(rowSums(tempK[,(origcolnumber+1):(dim(tempK)[2]-1)] == "1/0")))))
   tempK <- mutate(tempK, hom0 = rowSums(tempK[,(origcolnumber+1):(dim(tempK)[2])-2] == "0/0"))
   tempK <- filter(tempK,(!(((hom1+het)==0))|((hom0+het)==0)))
   write.table(headerrows,(paste(basename,".",parameters[2,1],"_",parameters[3,1],".",parameters[4,1],"_",parameters[6,1],".HWE.",popnames[k],".pop.vcf",sep="")),quote=FALSE,row.names=FALSE,col.names=FALSE)  
   write_delim(tempK[,1:origcolnumber],(paste(basename,".",parameters[2,1],"_",parameters[3,1],".",parameters[4,1],"_",parameters[6,1],".HWE.",popnames[k],".pop.vcf",sep="")),delim="\t",append=TRUE,col_names=TRUE)    
}
