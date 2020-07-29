#GBS_SNP_filter v1.17
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
temp <- temp %>% mutate_at(vars(10:dim(temp)[2]), .funs = list(cov = ~gsub("./.","0",gsub(":.*","",gsub("^.*?:","", . ))))) 
temp <- temp %>%  mutate_at(vars((origcolnumber+1):(dim(temp)[2])),list(~as.numeric(.)))

popmap <- read.table("popmap.txt",header=FALSE,stringsAsFactors=FALSE)
popnames <- unique(popmap[,2])

SNP_record <- NULL

removepops <- NULL

for (k in 1:length(popnames)) {
  if ((paste(basename,".",parameters[2,1],"_",parameters[3,1],".",parameters[4,1],"_",parameters[6,1],".HWE.",popnames[k],".pop.ld",sep="")) %in% list.files()) { #1A: Does this file exist (e.g. have all the sample been ditched due to SNP filters
    tempk <- read_table2((paste(basename,".",parameters[2,1],"_",parameters[3,1],".",parameters[4,1],"_",parameters[6,1],".HWE.",popnames[k],".pop.ld",sep="")),col_names=TRUE)[,1:7]
    tempk <- select(tempk,c(SNP_A,SNP_B,R2))
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

if(all(temp$ID==".")) {
  missing_ID <- TRUE
} else {
  missing_ID <- FALSE
}

while (j <= SNP_length) {
  # Getting the position of the SNPs in LD in the vcf
  if(missing_ID) {
    vcfpos1 <- which(temp$`#CHROM`==gsub(":.*","",SNP_record$Locus_1[j]) & temp$`POS`==gsub(".*:","",SNP_record$Locus_1[j]))
    vcfpos2 <- which(temp$`#CHROM`==gsub(":.*","",SNP_record$Locus_2[j]) & temp$`POS`==gsub(".*:","",SNP_record$Locus_2[j]))
  } else {
    vcfpos1 <- which(temp$ID==SNP_record$Locus_1[j])
    vcfpos2 <- which(temp$ID==SNP_record$Locus_2[j])
  }
  
  # Counting number of samples without data for each SNP (will ditch the SNP with more missing samples)
  zero_one_count <- length(which(temp[vcfpos1,(origcolnumber+1):(dim(temp)[2])]==0))
  zero_two_count <- length(which(temp[vcfpos2,(origcolnumber+1):(dim(temp)[2])]==0))
  
  # If Locus_1 SNP has more missing data than Locus_2
  if (zero_one_count > zero_two_count) {
    # Remove vcfpos1 from the vcf file
    temp <- temp[-vcfpos1,]
    # Record this snp as the one removed
    SNP_record[j,(dim(SNP_record)[2])] <- SNP_record[j,1]
    # Find the rows in SNP_record which also have the removed SNP
    todelete <- c((which(SNP_record[,1] %in% SNP_record[j,1])),(which(SNP_record[,2] %in% SNP_record[j,1])))     # Making sure that there would be something left in SNP_record if we deleted these rows (otherwise no point - we are at the end of the loop)
    if(length(todelete[(!(todelete <= j))])>0) {
        # Modifying todelete to not delete the row we've just recorded the LD SNP removed from
        todelete <- todelete[(!(todelete <= j))]
        # Removing these rows from SNP_record (b/c we've already deleted one of the SNPs in LD)
        SNP_record <- SNP_record[-todelete,]
    }
    # Moving on to the next row
    j <- j + 1  
  } else {
    # If Locus_2 SNP has more missing data than Locus_2
    if (zero_two_count > zero_one_count) {
      # Remove vcfpos2 from the vcf file
      temp <- temp[-vcfpos2,]
      # Record this snp as the one removed
      SNP_record[j,(dim(SNP_record)[2])] <- SNP_record[j,2] 
      # Find the rows in SNP_record which also have the removed SNP
      todelete <- c((which(SNP_record[,1] %in% SNP_record[j,2])),(which(SNP_record[,2] %in% SNP_record[j,2])))    
      # Making sure that there would be something left in SNP_record if we deleted these rows (otherwise no point - we are at the end of the loop)
      if(length(todelete[(!(todelete <= j))])>0) {
          # Modifying todelete to not delete the row we've just recorded the LD SNP removed from
          todelete <- todelete[(!(todelete <= j))]
          # Removing these rows from SNP_record (b/c we've already deleted one of the SNPs in LD)
          SNP_record <- SNP_record[-todelete,]
      }
      # Moving on to the next row
      j <- j + 1  
    } else {
      # If both SNPs have equal number of missing samples, will go by average coverage instead
      if (zero_two_count==zero_one_count) {
        zero_one_count <- mean(as.numeric(temp[vcfpos1,(origcolnumber+1):(dim(temp)[2])]))
        zero_two_count <- mean(as.numeric(temp[vcfpos2,(origcolnumber+1):(dim(temp)[2])]))
         # if Locus_1 SNP has lower average coverage than Locus_2
         if (zero_one_count < zero_two_count) {
            # Remove vcfpos1 from the vcf file
             temp <- temp[-vcfpos1,]
             # Record this snp as the one removed
             SNP_record[j,(dim(SNP_record)[2])] <- SNP_record[j,1]
             # Find the rows in SNP_record which also have the removed SNP
             todelete <- c((which(SNP_record[,1] %in% SNP_record[j,1])),(which(SNP_record[,2] %in% SNP_record[j,1])))    
             # Making sure that there would be something left in SNP_record if we deleted these rows (otherwise no point - we are at the end of the loop)
             if(length(todelete[(!(todelete <= j))])>0) {
                # Modifying todelete to not delete the row we've just recorded the LD SNP removed from
                todelete <- todelete[(!(todelete <= j))]
                # Removing these rows from SNP_record (b/c we've already deleted one of the SNPs in LD)
                SNP_record <- SNP_record[-todelete,]
             }
             # Moving on to the next row
            j <- j + 1  
          } else {
            # if Locus_2 SNP has lower average coverage than Locus_1
            if (zero_two_count < zero_one_count) {
              # Remove vcfpos2 from the vcf file
              temp <- temp[-vcfpos2,]
              # Record this snp as the one removed
              SNP_record[j,(dim(SNP_record)[2])] <- SNP_record[j,2]
              # Find the rows in SNP_record which also have the removed SNP
              todelete <- c((which(SNP_record[,1] %in% SNP_record[j,2])),(which(SNP_record[,2] %in% SNP_record[j,2])))
              # Making sure that there would be something left in SNP_record if we deleted these rows (otherwise no point - we are at the end of the loop)
              if(length(todelete[(!(todelete <= j))])>0) {
                 # Modifying todelete to not delete the row we've just recorded the LD SNP removed from
                 todelete <- todelete[(!(todelete <= j))]
                 # Removing these rows from SNP_record (b/c we've already deleted one of the SNPs in LD)
                 SNP_record <- SNP_record[-todelete,]
              }
              # Moving on to the next row
              j <- j + 1  
            } else {
              # If the SNPs are equal in coverage
              if (zero_two_count==zero_one_count) {
                  # Arbritrarily chosing to remove the second SNP
                  temp <- temp[-vcfpos2,]
                  # Record this snp as the one removed
                  SNP_record[j,(dim(SNP_record)[2])] <- SNP_record[j,2]
                  # Find the rows in SNP_record which also have the removed SNP
                  todelete <- c((which(SNP_record[,1] %in% SNP_record[j,2])),(which(SNP_record[,2] %in% SNP_record[j,2])))
                  # Making sure that there would be something left in SNP_record if we deleted these rows (otherwise no point - we are at the end of the loop)
                  if(length(todelete[(!(todelete <= j))])>0) {
                    # Modifying todelete to not delete the row we've just recorded the LD SNP removed from
                    todelete <- todelete[(!(todelete <= j))]
                    # Removing these rows from SNP_record (b/c we've already deleted one of the SNPs in LD)
                    SNP_record <- SNP_record[-todelete,]
                  }
                  # Moving on to the next row
                  j <- j + 1              
               }
             }
          }
       }  
    }
  }
  print(paste("Up to ",j," out of ",dim(SNP_record)[1]," pairwise LD comparisons",sep=""))
  SNP_length <- dim(SNP_record)[1]
  if (SNP_length==1) {
    break
  }
}    

write(format(Sys.time(),usetz = TRUE),logfilename,append=TRUE)
write(paste(dim(SNP_record)[1]," loci will be removed as they were in linkage with another locus in at least ",parameters[6,1]," populations at an Rsq of >=",parameters[5,1],sep=""),logfilename,append=TRUE)
write(paste("These loci will be listed in ",basename,".",parameters[2,1],"_",parameters[3,1],".",parameters[4,1],"_",parameters[6,1],".HWE.",parameters[5,1],".rsq",sep=""),logfilename,append=TRUE)
write.table(SNP_record,paste(basename,".",parameters[2,1],"_",parameters[3,1],".",parameters[4,1],"_",parameters[6,1],".HWE.",parameters[5,1],".rsq",sep=""),row.names=FALSE,col.names=TRUE,quote=FALSE)  
write.table(headerrows,(paste(basename,".",parameters[2,1],"_",parameters[3,1],".",parameters[4,1],"_",parameters[6,1],".HWE.",parameters[5,1],".ld.vcf",sep="")),quote=FALSE,row.names=FALSE,col.names=FALSE)
write_delim(temp[,1:origcolnumber],(paste(basename,".",parameters[2,1],"_",parameters[3,1],".",parameters[4,1],"_",parameters[6,1],".HWE.",parameters[5,1],".ld.vcf",sep="")),delim="\t",append=TRUE,col_names=TRUE)    
write(format(Sys.time(),usetz = TRUE),logfilename,append=TRUE)  
write(paste("Following this filtering ",basename,".",parameters[2,1],"_",parameters[3,1],".",parameters[4,1],"_",parameters[6,1],".HWE.",parameters[5,1],".ld.vcf has been written out, containing ",(dim(temp)[1])," SNPs and ", (origcolnumber-9), " samples",sep=""),logfilename,append=TRUE)   

print("GBS_SNP_filter_rsq.R has the following warnings():")
warnings()
print("If no warnings printed, none thrown by GBS_SNP_filter_rsq.R")
