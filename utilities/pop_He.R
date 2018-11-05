pop_He <- function(vcf_name,min_sample_size) {

  if (!require('vcfR')) install.packages('vcfR'); library('vcfR')
  vcf <- read.vcfR(vcf_name, verbose = FALSE )
  popsum <- gt.to.popsum(vcf)
  basename <- gsub(".vcf","",vcf_name)
  popsum <- popsum[(which(popsum[,1]>=min_sample_size)),]
  print(paste("The average He across all loci in ",vcf_name," is: ",mean(popsum[,3]),sep=""))
  print(paste("The effective number of alleles averaged across all loci in ",vcf_name," is: ",mean(popsum[,4]),sep=""))
  print("The per locus values have been written out to:")
  print(paste(basename,".",min_sample_size,".He.txt",sep=""))
  write.table(popsum,(paste(basename,".",min_sample_size,".He.txt",sep="")),quote=FALSE,row.names=FALSE)  

}  
