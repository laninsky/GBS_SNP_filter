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
temp <- temp %>% mutate_at(vars(10:dim(temp)[2]), .funs = funs(cov = gsub(":.*","",gsub("^.*?:","", . )))) 
temp <- temp %>%  mutate_at(vars((origcolnumber+1):(dim(temp)[2])),funs(as.numeric))

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

SNP_record <- filter(SNP_record, (dim(SNP_record)[2]-2)-rowSums(is.na(SNP_record[,3:(dim(SNP_record)[2])]))>=as.numeric(parameters[6,1]))

SNP_record <- cbind(SNP_record, (matrix(NA,ncol=1,nrow=dim(SNP_record)[1])))
names(SNP_record)[(dim(SNP_record)[2])] <- "removed"

j <- 1
SNP_length <- dim(SNP_record)[1]

while (j <= SNP_length) {
  zero_one_count <- sum(temp[(which(temp$`#CHROM` %in% SNP_record[j,1])),((origcolnumber+1):(dim(temp)[2]))]==0)
  zero_two_count <- sum(temp[(which(temp$`#CHROM` %in% SNP_record[j,2])),((origcolnumber+1):(dim(temp)[2]))]==0)
  if (zero_one_count > zero_two_count) {
    temp <- temp[-(which(temp$`#CHROM` %in% SNP_record[j,1])),]
    SNP_record[j,(dim(SNP_record)[2])] <- SNP_record[j,1]    
    todelete <- c((which(SNP_record[,1] %in% SNP_record[j,1])),(which(SNP_record[,2] %in% SNP_record[j,1])))    
    if(length(todelete[(!(todelete <= j))])>0) {
        todelete <- todelete[(!(todelete <= j))]
        SNP_record <- SNP_record[-todelete,]
    }
    j <- j + 1  
  } else {
    if (zero_two_count > zero_one_count) {
      temp <- temp[-(which(temp$`#CHROM` %in% SNP_record[j,2])),]
      SNP_record[j,(dim(SNP_record)[2])] <- SNP_record[j,2]    
      todelete <- c((which(SNP_record[,1] %in% SNP_record[j,2])),(which(SNP_record[,2] %in% SNP_record[j,2])))    
      if(length(todelete[(!(todelete <= j))])>0) {
          todelete <- todelete[(!(todelete <= j))]
          SNP_record <- SNP_record[-todelete,]
      }
      j <- j + 1  
    } else {
      if (zero_two_count==zero_one_count) {
        zero_one_count <- sum(temp[(which(temp$`#CHROM` %in% SNP_record[j,1])),((origcolnumber+1):(dim(temp)[2]))])
        zero_two_count <- sum(temp[(which(temp$`#CHROM` %in% SNP_record[j,2])),((origcolnumber+1):(dim(temp)[2]))])         
         if (zero_one_count < zero_two_count) {
            temp <- temp[-(which(temp$`#CHROM` %in% SNP_record[j,1])),]
            SNP_record[j,(dim(SNP_record)[2])] <- SNP_record[j,1]    
            todelete <- c((which(SNP_record[,1] %in% SNP_record[j,1])),(which(SNP_record[,2] %in% SNP_record[j,1])))    
            if(length(todelete[(!(todelete <= j))])>0) {
              todelete <- todelete[(!(todelete <= j))]
              SNP_record <- SNP_record[-todelete,]
            }
            j <- j + 1  
          } else {
            if (zero_two_count < zero_one_count) {
              temp <- temp[-(which(temp$`#CHROM` %in% SNP_record[j,2])),]
              SNP_record[j,(dim(SNP_record)[2])] <- SNP_record[j,2]    
              todelete <- c((which(SNP_record[,1] %in% SNP_record[j,2])),(which(SNP_record[,2] %in% SNP_record[j,2])))    
              if(length(todelete[(!(todelete <= j))])>0) {
                 todelete <- todelete[(!(todelete <= j))]
                 SNP_record <- SNP_record[-todelete,]
              }
              j <- j + 1  
            } else {
              if (zero_two_count==zero_one_count) {
                 temp <- temp[-(which(temp$`#CHROM` %in% SNP_record[j,2])),]
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
  #OK, SO MAYBE WE WRITE OUT THIS TABLE. DOUBLE CHECK THINGS ARE NOT WEIRD, BUT ACTUALLY THAT SHOULD BE IT?
  # CHECK NOTHING IS IN THE LEFT AND RIGHT COLUMN
}    


#Probably want to write out this SNP list
#Then want to use it to filter out the SNPs that are in linkage, keeping the SNP with the highest coverage. When the SNP is ditched, need to check SNP_record and get rid of any rows where the ditched SNP also occurs.