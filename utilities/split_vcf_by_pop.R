numberofheaders <- 10
vcf_name <- "hihi_prev_filter.0.85_0.9.0.05_3.HWE.0.5.ld.vcf"

split_vcf_by_ppop <- function(vcf_name,numberofheaders) {
  if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
  if (!require('readr')) install.packages('readr'); library('readr')
  if (!require('stringr')) install.packages('stringr'); library('stringr')
  temp <- read_tsv(vcf_name,col_names=TRUE,skip=numberofheaders)
  headerrows <- readLines(vcf_name,numberofheaders)

  popmap <- read.table("popmap.txt",header=FALSE,stringsAsFactors=FALSE)
  popnames <- unique(popmap[,2])
  for (k in 1:length(popnames)) {
    tempK <- select(temp, c((1:9),which(names(temp) %in% (popmap[(which(popmap[,2]==popnames[k])),1]))))
    origcolnumber <- dim(tempK)[2]
    tempK <- mutate_at(tempK,vars(10:origcolnumber),.funs = funs(genotype = gsub(":.*","", . )))     
    tempK <- mutate(tempK, hom1 = rowSums(tempK[,(origcolnumber+1):(dim(tempK)[2])] == "1/1"))
    tempK <- mutate(tempK, het = ((rowSums(tempK[,(origcolnumber+1):(dim(tempK)[2]-1)] == "0/1")+(rowSums(tempK[,(origcolnumber+1):(dim(tempK)[2]-1)] == "1/0")))))
    tempK <- mutate(tempK, hom0 = rowSums(tempK[,(origcolnumber+1):(dim(tempK)[2])-2] == "0/0"))
    tempK <- filter(tempK,(!(((hom1+het)==0))|((hom0+het)==0)))
    write.table(headerrows,(paste(basename,".",popnames[k],".pop.vcf",sep="")),quote=FALSE,row.names=FALSE,col.names=FALSE)  
    write_delim(tempK[,1:origcolnumber],paste(basename,".",popnames[k],".pop.vcf",sep="")),delim="\t",append=TRUE,col_names=TRUE)    
  }
}  
