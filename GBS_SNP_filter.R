library(tidyverse)
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
}  

if (!((paste(basename,".",parameters[2,1],"_",parameters[3,1],".HWE.vcf",sep="")) %in% filelist)) {
  if (!((paste(basename,".HWE",sep="")) %in% filelist)) { # If locus specific HWE values have not already been printed out  
    popmap <- read.table("popmap.txt",header=FALSE,stringsAsFactors=FALSE)
    popnames <- unique(popmap[,2])
    hwetable <- matrix(c("snp",popnames),nrow=1)
    write.table(hwetable,(paste(basename,".HWE",sep="")),quote=FALSE,row.names=FALSE,col.names=FALSE)  
    hwetablebin <- hwetable
    for (i in 1:(dim(temp)[1])) { #5A: for each SNP
      temprow <- matrix(c(temp[i,1],popnames),nrow=1)
      temptemp <- temp[i,1:origcolnumber]      
      for (k in 1:length(popnames)) { #8A: for each population
        temptemppop <- select(temptemp, which(names(temptemp) %in% (popmap[(which(popmap[,2]==popnames[k])),1])))
        temptemppop <- mutate_at(temptemppop,vars(1:dim(temptemppop)[2]),funs(gsub(":.*","", . )))
        tempmatrix <- matrix(0,ncol=2,nrow=3)
        tempmatrix[1,1] <- length(which(temptemppop[1,]=="0/0"))
        tempmatrix[2,1] <- length(which((temptemppop[1,]=="0/1" | temptemppop[1,]=="1/0")))
        tempmatrix[3,1] <- length(which(temptemppop[1,]=="1/1"))
        tempmatrix[1,2] <- ((((2*tempmatrix[1,1])+tempmatrix[2,1])/(2*sum(tempmatrix[,1])))^2)*sum(tempmatrix[,1])
        tempmatrix[3,2] <- ((((2*tempmatrix[3,1])+tempmatrix[2,1])/(2*sum(tempmatrix[,1])))^2)*sum(tempmatrix[,1])
        tempmatrix[2,2] <- 2*(((2*tempmatrix[1,1])+tempmatrix[2,1])/(2*sum(tempmatrix[,1])))*(((2*tempmatrix[3,1])+tempmatrix[2,1])/(2*sum(tempmatrix[,1])))*sum(tempmatrix[,1])
        if (sum(tempmatrix[,1])==0) {
          temprow[1,(k+1)] <- "NaN"
        } else {  
          temprow[1,(k+1)] <- suppressWarnings(fisher.test(tempmatrix)$p.value)
        }  
      } #8B
      if(sum(temprow[2:length(temprow)]<as.numeric(parameters[4,1]))>as.numeric(parameters[6,1])) {
        hwetablebin <- rbind(hwetablebin,temprow)
      }
      write.table(temprow,(paste(basename,".HWE",sep="")),quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
      print(paste("Up to ",i," out of ",(dim(temp)[1]), " loci",sep=""))  
    }     
    write(format(Sys.time(),usetz = TRUE),logfilename,append=TRUE)
    write(paste("The following loci (p-values given) will be removed as more than ",parameters[6,1]," populations had a HWE p-value of <",parameters[4,1],sep=""),logfilename,append=TRUE)
    write.table(hwetablebin,logfilename,append=TRUE,row.names=FALSE,col.names=FALSE)  
    temp <- temp %>% filter(., (!(`#CHROM` %in% (hwetablebin[2:(dim(hwetablebin)[2]),1]))))
    write.table(headerrows,(paste(basename,".",parameters[2,1],"_",parameters[3,1],".HWE.vcf",sep="")),quote=FALSE,row.names=FALSE,col.names=FALSE)
    write_delim(temp[,1:origcolnumber],(paste(basename,".",parameters[2,1],"_",parameters[3,1],".HWE.vcf",sep="")),delim="\t",append=TRUE,col_names=TRUE)    
    write(format(Sys.time(),usetz = TRUE),logfilename,append=TRUE)  
    write(paste("Following this filtering ",basename,".",parameters[2,1],"_",parameters[3,1],".HWE.vcf has been written out, containing ",(dim(temp)[1])," SNPs and ", (origcolnumber-9), " samples",sep=""),logfilename,append=TRUE)  
  } else { #what to do if *.HWE does exist and you can use it to filter things    
    hwetable <- read.table((paste(basename,".HWE",sep="")),header=TRUE,stringsAsFactors=FALSE)
    hwetablebin <- NULL
    for (i in 1:(dim(hwetable)[1])) {
        if(sum(hwetable[i,2:(dim(hwetable)[2])]<as.numeric(parameters[4,1]))>as.numeric(parameters[6,1])) {
            hwetablebin <- rbind(hwetablebin,hwetable[i,])
        }  
    }
    write(format(Sys.time(),usetz = TRUE),logfilename,append=TRUE)
    write(paste("The following loci (p-values given) will be removed as more than ",parameters[6,1]," populations had a HWE p-value of <",parameters[4,1],sep=""),logfilename,append=TRUE)
    write.table(hwetablebin,logfilename,append=TRUE,row.names=FALSE)
    temp <- temp %>% filter(., (!(`#CHROM` %in% (hwetablebin[1:(dim(hwetablebin)[2]),1]))))
    write.table(headerrows,(paste(basename,".",parameters[2,1],"_",parameters[3,1],".HWE.vcf",sep="")),quote=FALSE,row.names=FALSE,col.names=FALSE)  
    write_delim(temp[,1:origcolnumber],(paste(basename,".",parameters[2,1],"_",parameters[3,1],".HWE.vcf",sep="")),delim="\t",append=TRUE,col_names=TRUE)    
    write(format(Sys.time(),usetz = TRUE),logfilename,append=TRUE)  
    write(paste("Following this filtering ",basename,".",parameters[2,1],"_",parameters[3,1],".HWE.vcf has been written out, containing ",(dim(temp)[1])," SNPs and ", (origcolnumber-9), " samples",sep=""),logfilename,append=TRUE)     
  }
} else {  
  headerrows <- read_tsv("header_row.txt",col_names=FALSE)
  numberofheaders <- dim(headerrows)[1]
  temp <- read_tsv((paste(basename,".",parameters[2,1],"_",parameters[3,1],".HWE.vcf",sep="")),col_names=TRUE,skip=numberofheaders)
  origcolnumber <- dim(temp)[2]
  temp <- temp %>% mutate_at(vars(10:dim(temp)[2]), .funs = funs(cov = gsub(":.*","",gsub("^.*?:","", . )))) 
  temp <- temp %>% mutate_at(vars((origcolnumber+1):(dim(temp)[2])),funs(as.numeric))
}

#RSQ is too computationally costly to do on the "full dataset". Instead, bring the HWE calculations up here,
# After filtering on this, then can do Rsq at the end. NEED TO ADD IN OPTIONS IF RSQ EXISTS

if (!((paste(basename,".",parameters[2,1],"_",parameters[3,1],"_",parameters[5,1],".rsq",sep="")) %in% filelist)) { #4A: if *.rsq doesn't exist, creating this
  otherrsq <- strsplit(list.files(pattern=".rsq"),".rsq")
  rsqfile <- NULL
  for (i in 1:length(otherrsq)) {
    eachfile <- unlist(strsplit(otherrsq[[i]],"_"))
    suppressWarnings( if(is.na(as.numeric(eachfile[length(eachfile)]))) {
      break
    })  
    if(as.numeric(eachfile[length(eachfile)])<=parameters[5,1]) {
      suppressWarnings( if(is.na(as.numeric(eachfile[length(eachfile)]))) {
        break
      })  
      if(as.numeric(eachfile[length(eachfile)-1])>=parameters[3,1]) {
        suppressWarnings( if(is.na(as.numeric(paste(unlist(strsplit(eachfile[length(eachfile)-2],".",fixed=TRUE))[2:3],collapse=".")))) {
          break
        })  
        if(as.numeric(paste(unlist(strsplit(eachfile[length(eachfile)-2],".",fixed=TRUE))[2:3],collapse="."))<=parameters[2,1]) {
          rsqfile <- paste(otherrsq[[i]],".rsq",sep="")
          break
        }
      }
    }
  }  
  if(!(is.null(rsqfile))) {  #5A If there are no other *rsq files that can be used...     
    popmap <- read.table("popmap.txt",header=FALSE,stringsAsFactors=FALSE)
    popnames <- unique(popmap[,2])
    LDbin <- matrix(c("snp1","snp2","popname","rsq"),nrow=1)
    write.table(LDbin,(paste(basename,".",parameters[2,1],"_",parameters[3,1],"_",parameters[5,1],".rsq",sep="")),quote=FALSE,row.names=FALSE,col.names=FALSE)    
    for (k in 1:length(popnames)) { #6A for each pop
      tempK <- select(temp, c(`#CHROM`,which(names(temp) %in% (popmap[(which(popmap[,2]==popnames[k])),1]))))
      tempK <- mutate_at(tempK,vars(2:dim(tempK)[2]),funs(gsub(":.*","", . )))     
      tempK <- mutate(tempK, hom1 = rowSums(tempK[,2:(dim(tempK)[2])] == "1/1"))
      tempK <- mutate(tempK, het = ((rowSums(tempK[,2:(dim(tempK)[2]-1)] == "0/1")+(rowSums(tempK[,2:(dim(tempK)[2]-1)] == "1/0")))))
      tempK <- mutate(tempK, hom0 = rowSums(tempK[,2:(dim(tempK)[2])-2] == "0/0"))
      tempK <- filter(tempK,(!(((hom1+het)==0))|((hom0+het)==0)))
      totalcomps <- (dim(tempK)[1]*(dim(tempK)[1]-1))/2
      culmcompls <- 0
      for (i in 1:(dim(tempK)[1]-1) { #7A for each locus
        zerocounts <- unlist(lapply(((i+1):(dim(tempK)[1])),function(x){
          tempmatrix <- matrix(0,ncol=3,nrow=3)
          tempmatrix[1,1] <- length(which(tempK[i,]=="0/0" & tempK[x,]=="0/0"))
          tempmatrix[2,1] <- length(which((tempK[i,]=="0/1" | tempK[i,]=="1/0") & tempK[x,]=="0/0"))
          tempmatrix[3,1] <- length(which(tempK[i,]=="1/1" & tempK[x,]=="0/0"))
          tempmatrix[1,2] <- length(which(tempK[i,]=="0/0" & (tempK[x,]=="0/1" | tempK[x,]=="1/0")))
          tempmatrix[2,2] <- length(which((tempK[i,]=="0/1" | tempK[i,]=="1/0") & (tempK[x,]=="0/1" | tempK[x,]=="1/0")))
          tempmatrix[3,2] <- length(which(tempK[i,]=="1/1" & (tempK[x,]=="0/1" | tempK[x,]=="1/0")))
          tempmatrix[1,3] <- length(which(tempK[i,]=="0/0" & tempK[x,]=="1/1"))
          tempmatrix[2,3] <- length(which((tempK[i,]=="0/1" | tempK[i,]=="1/0") & tempK[x,]=="1/1"))
          tempmatrix[3,3] <- length(which(tempK[i,]=="1/1" & tempK[x,]=="1/1"))
          twobytwo <- matrix(0,nrow=2,ncol=2)
          twobytwo[1,1] <- 2*tempmatrix[1,1]+tempmatrix[2,1]+tempmatrix[1,2]
          twobytwo[2,1] <- 2*tempmatrix[3,1]+tempmatrix[3,2]+tempmatrix[2,1]
          twobytwo[1,2] <- 2*tempmatrix[1,3]+tempmatrix[1,2]+tempmatrix[2,3]
          twobytwo[2,2] <- 2*tempmatrix[3,3]+tempmatrix[3,2]+tempmatrix[2,3]
          oddsratio <- (twobytwo[1,1]/twobytwo[2,1])/(twobytwo[1,2]/twobytwo[2,2])
          if(is.na(oddsratio)) {
            oddsratio <- 0
          }
          if(!(is.infinite(oddsratio))) {
            twobytwo[1,1] <- twobytwo[1,1] + tempmatrix[2,2]*oddsratio/(1+oddsratio)
            twobytwo[2,1] <- twobytwo[2,1] + tempmatrix[2,2]*1/(1+oddsratio)  
            twobytwo[1,2] <- twobytwo[1,2] + tempmatrix[2,2]*1/(1+oddsratio) 
            twobytwo[2,2] <- twobytwo[2,2] + tempmatrix[2,2]*oddsratio/(1+oddsratio)
          }  
          return(((twobytwo[1,1]*twobytwo[2,2]-twobytwo[1,2]*twobytwo[2,1])^2)/(sum(twobytwo[,1])*sum(twobytwo[,2])*sum(twobytwo[1,])*sum(twobytwo[2,])))
        }))
        zerocountpos <- which(zerocounts>as.numeric(parameters[5,1]))+i
        LDbintemp <- matrix(NA,nrow=length(zerocountpos),ncol=4)
        LDbintemp[,1] <- as.matrix(tempK[i,1])
        LDbintemp[,2] <- as.matrix(tempK[zerocountpos,1])
        LDbintemp[,3] <- popnames[k]
        LDbintemp[,4] <- zerocounts[zerocountpos-i]       
        write.table(LDbintemp,(paste(basename,".",parameters[2,1],"_",parameters[3,1],"_",parameters[5,1],".rsq",sep="")),quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
        culmcompls <- culmcompls+((dim(tempK)[1]-i))
        paste((culmcompls/totalcomps*100),"% through ",totalcomps, " pairwise SNP comparisons for ", popnames[k]," (",k," out of ",length(popnames)," populations)",sep="")
        
      } #7B for each SNP
    } #8B for each population       
  } else { #5AB
       #what to do if there is an rsq table you can read in
  } #5B   
} else { #4AB: if the exact file needed exists
    #what to do if there is the exact rsq table you need that can read in
} #4B  
    
    
## TOO SLOW. FIGURE OUT HOW TO VECTORIZE
if (!((paste(basename,".rsq",sep="")) %in% filelist)) { #4A: if *.rsq doesn't exist, creating this
  popmap <- read.table("popmap.txt",header=FALSE,stringsAsFactors=FALSE)
  popnames <- unique(popmap[,2])
  tablerow <- matrix(c("snp1","snp2",popnames),nrow=1)
  ldbin <- tablerow
  write.table(tablerow,(paste(basename,".rsq",sep="")),quote=FALSE,row.names=FALSE,col.names=FALSE) 
  for (i in 1:(dim(temp)[1]-1)) { #5A: for SNP1
    for (j in (i+1):(dim(temp)[1])) { #6A for SNP2
      temprow <- matrix(c(temp[i,1],temp[j,1],popnames),nrow=1)
      temptemp <- temp[c(i,j),1:origcolnumber]      
      for (k in 1:length(popnames)) { #7A: for each population

        } #7B
      write.table(temprow,(paste(basename,".rsq",sep="")),quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
      if(sum(unlist(temprow[1,3:length(temprow)])>as.numeric(parameters[5,1]),na.rm=TRUE)>as.numeric(parameters[6,1])) {
        ldbin <- rbind(ldbin,temprow)
      }  
    } #6B
  } #5B  
} #4B  
 
  
