popmap <- read.table("popmap.txt",header=FALSE,stringsAsFactors=FALSE)
popnames <- unique(popmap[,2])
for (k in 1:length(popnames)) {
  # Probably need a fix here too
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
