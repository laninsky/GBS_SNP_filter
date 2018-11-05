pop_He <- function(vcf_name,min_sample_size) {

  if (!require('vcfR')) install.packages('vcfR'); library('vcfR')
  vcf <- read.vcfR(vcf_name, verbose = FALSE )
  popsum <- gt.to.popsum(vcf)
  basename <- gsub(".vcf","",vcf_name)
  popsum <- popsum[(which(popsum[,1]>=min_sample_size)),]
  write.table(popsum,(paste(basename,".",popnames[k],".",min_sample_size,".He.txt",sep="")),quote=FALSE,row.names=FALSE)  

}  
