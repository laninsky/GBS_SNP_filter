# Loading in relevant packages
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('readr')) install.packages('readr'); library('readr')
if (!require('stringr')) install.packages('stringr'); library('stringr')
if (!require('tibble')) install.packages('tibble'); library('tibble')
#library(tidyverse) # I've had issues loading library(tidyverse) and R crashing using the sbatch syste
parameters <- read.table("GBS_SNP_filter.txt",header=FALSE,stringsAsFactors=FALSE,comment.char = "")
basename <- gsub(".vcf","",parameters[1,1])
filelist <- list.files()
logfilename <- paste(basename,".log",sep="")

headerrows <- read_tsv("header_row.txt",col_names=FALSE)
numberofheaders <- dim(headerrows)[1]

temp <- read_tsv((paste(basename,".",parameters[2,1],"_",parameters[3,1],".",parameters[4,1],"_",parameters[6,1],".HWE.vcf",sep="")),col_names=TRUE,skip=numberofheaders)
origcolnumber <- dim(temp)[2]
temp <- temp %>% mutate_at(vars(10:dim(temp)[2]), .funs = funs(cov = gsub("./.","0",gsub(":.*","",gsub("^.*?:","", . ))))) 
temp <- temp %>%  mutate_at(vars((origcolnumber+1):(dim(temp)[2])),funs(as.numeric))

popmap <- read.table("popmap.txt",header=FALSE,stringsAsFactors=FALSE)
popnames <- unique(popmap[,2])

SNP_record <- NULL

removepops <- NULL

for (k in 1:length(popnames)) {
  if ((paste(basename,".",parameters[2,1],"_",parameters[3,1],".",parameters[4,1],"_",parameters[6,1],".HWE.",popnames[k],".pop.ld",sep="")) %in% list.files()) { #1A: Does this file exist (e.g. have all the sample been ditched due to SNP filters
    tempk <- read_table2((paste(basename,".",parameters[2,1],"_",parameters[3,1],".",parameters[4,1],"_",parameters[6,1],".HWE.",popnames[k],".pop.ld",sep="")),col_names=TRUE)[,1:7]
    tempk <- select(tempk,c(SNP_A,SNP_B,R2))
    tempk <- mutate_at(tempk,vars(SNP_A,SNP_B),funs(gsub(":.*","", . )))
    if(is.null(SNP_record)) {
      SNP_record <- tempk    
    }  else {
      SNP_record <- full_join(SNP_record,tempk,by=c("SNP_A","SNP_B"))
    }
  } else { #1AB: What to do if that population does have no samples
    write(format(Sys.time(),usetz = TRUE),logfilename,append=TRUE)
    write(paste(popnames[k]," does not have any samples remaining, potentially because they were removed due to having too many missing SNPs",sep=""),logfilename,append=TRUE)
    removepops <- c(removepops, k)
  } #1B  
}

# Modified because of https://github.com/laninsky/GBS_SNP_filter/issues/1 - I was assuming populations would be dropped!
if(length(removepops)>0) {
  names(SNP_record) <- c("Locus_1","Locus_2",popnames[-removepops])
} else {
  names(SNP_record) <- c("Locus_1","Locus_2",popnames)
}

SNP_record <- filter(SNP_record, (dim(SNP_record)[2]-2)-rowSums(is.na(SNP_record[,3:(dim(SNP_record)[2])]))>=as.numeric(parameters[6,1]))

SNP_record <- cbind(SNP_record, (matrix(NA,ncol=1,nrow=dim(SNP_record)[1])))
names(SNP_record)[(dim(SNP_record)[2])] <- "removed"

j <- 1
SNP_length <- dim(SNP_record)[1]
locusname <- temp %>% select(!!parameters[7,1])
temp <- add_column(temp,as.matrix(locusname)[,1],.before=TRUE)
names(temp)[1] <- "locusname"  

while (j <= SNP_length) {
  zero_one_count <- sum(temp[(which(temp$locusname %in% SNP_record[j,1])),((origcolnumber+2):(dim(temp)[2]))]==0)
  zero_two_count <- sum(temp[(which(temp$locusname %in% SNP_record[j,2])),((origcolnumber+2):(dim(temp)[2]))]==0)
  
  if (zero_one_count > zero_two_count) {
    temp <- temp[-(which(temp$locusname %in% SNP_record[j,1])),]
    SNP_record[j,(dim(SNP_record)[2])] <- SNP_record[j,1]    
    todelete <- c((which(SNP_record[,1] %in% SNP_record[j,1])),(which(SNP_record[,2] %in% SNP_record[j,1])))    
    if(length(todelete[(!(todelete <= j))])>0) {
        todelete <- todelete[(!(todelete <= j))]
        SNP_record <- SNP_record[-todelete,]
    }
    j <- j + 1  
  } else {
    if (zero_two_count > zero_one_count) {
      temp <- temp[-(which(temp$locusname %in% SNP_record[j,2])),]
      SNP_record[j,(dim(SNP_record)[2])] <- SNP_record[j,2]    
      todelete <- c((which(SNP_record[,1] %in% SNP_record[j,2])),(which(SNP_record[,2] %in% SNP_record[j,2])))    
      if(length(todelete[(!(todelete <= j))])>0) {
          todelete <- todelete[(!(todelete <= j))]
          SNP_record <- SNP_record[-todelete,]
      }
      j <- j + 1  
    } else {
      if (zero_two_count==zero_one_count) {
        zero_one_count <- sum(temp[(which(temp$locusname %in% SNP_record[j,1])),((origcolnumber+2):(dim(temp)[2]))])
        zero_two_count <- sum(temp[(which(temp$locusname %in% SNP_record[j,2])),((origcolnumber+2):(dim(temp)[2]))])         
         if (zero_one_count < zero_two_count) {
            temp <- temp[-(which(temp$locusname %in% SNP_record[j,1])),]
            SNP_record[j,(dim(SNP_record)[2])] <- SNP_record[j,1]    
            todelete <- c((which(SNP_record[,1] %in% SNP_record[j,1])),(which(SNP_record[,2] %in% SNP_record[j,1])))    
            if(length(todelete[(!(todelete <= j))])>0) {
              todelete <- todelete[(!(todelete <= j))]
              SNP_record <- SNP_record[-todelete,]
            }
            j <- j + 1  
          } else {
            if (zero_two_count < zero_one_count) {
              temp <- temp[-(which(temp$locusname %in% SNP_record[j,2])),]
              SNP_record[j,(dim(SNP_record)[2])] <- SNP_record[j,2]    
              todelete <- c((which(SNP_record[,1] %in% SNP_record[j,2])),(which(SNP_record[,2] %in% SNP_record[j,2])))    
              if(length(todelete[(!(todelete <= j))])>0) {
                 todelete <- todelete[(!(todelete <= j))]
                 SNP_record <- SNP_record[-todelete,]
              }
              j <- j + 1  
            } else {
              if (zero_two_count==zero_one_count) {
                 temp <- temp[-(which(temp$locusname %in% SNP_record[j,2])),]
                 SNP_record[j,(dim(SNP_record)[2])] <- SNP_record[j,2]    
                 todelete <- c((which(SNP_record[,1] %in% SNP_record[j,2])),(which(SNP_record[,2] %in% SNP_record[j,2])))    
                 if(length(todelete[(!(todelete <= j))])>0) {
                   todelete <- todelete[(!(todelete <= j))]
                   SNP_record <- SNP_record[-todelete,]
                 }
                 j <- j + 1              
               }
             }
          }
       }  
    }
  }
  print(paste("Up to ",j," out of ",dim(SNP_record)[1]," pairwise LD comparisons",sep=""))
  SNP_length <- dim(SNP_record)[1]
  if(is.na(SNP_record[(j-1),8])) {
    break
  }
}    

write(format(Sys.time(),usetz = TRUE),logfilename,append=TRUE)
write(paste(dim(SNP_record)[1]," loci will be removed as they were in linkage with another locus in at least ",parameters[6,1]," populations at an Rsq of >=",parameters[5,1],sep=""),logfilename,append=TRUE)
write(paste("These loci will be listed in ",basename,".",parameters[2,1],"_",parameters[3,1],".",parameters[4,1],"_",parameters[6,1],".HWE.",parameters[5,1],".rsq",sep=""),logfilename,append=TRUE)
write.table(SNP_record,paste(basename,".",parameters[2,1],"_",parameters[3,1],".",parameters[4,1],"_",parameters[6,1],".HWE.",parameters[5,1],".rsq",sep=""),row.names=FALSE,col.names=TRUE,quote=FALSE)  
write.table(headerrows,(paste(basename,".",parameters[2,1],"_",parameters[3,1],".",parameters[4,1],"_",parameters[6,1],".HWE.",parameters[5,1],".ld.vcf",sep="")),quote=FALSE,row.names=FALSE,col.names=FALSE)
write_delim(temp[,2:(origcolnumber+1)],(paste(basename,".",parameters[2,1],"_",parameters[3,1],".",parameters[4,1],"_",parameters[6,1],".HWE.",parameters[5,1],".ld.vcf",sep="")),delim="\t",append=TRUE,col_names=TRUE)    
write(format(Sys.time(),usetz = TRUE),logfilename,append=TRUE)  
write(paste("Following this filtering ",basename,".",parameters[2,1],"_",parameters[3,1],".",parameters[4,1],"_",parameters[6,1],".HWE.",parameters[5,1],".ld.vcf has been written out, containing ",(dim(temp)[1])," SNPs and ", (origcolnumber-9), " samples",sep=""),logfilename,append=TRUE)   
