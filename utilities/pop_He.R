pop_He <- function(vcf_name) {  
if (!require('vcfR')) install.packages('vcfR'); library('vcfR')
if (!require('reshape')) install.packages('reshape'); library('reshape')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')

vcf <- read.vcfR(vcf_name, verbose = FALSE )
gt <- extract.gt(vcf)
gt <- as.data.frame(t(gt), stringsAsFactors = FALSE)
gt[gt == "."] <- NA

snp_geno <- as.data.frame(sapply(gt,function(x){as.numeric(gsub("/","",x,fixed=TRUE))}))
row.names(snp_geno) <- row.names(gt)

popmap <- read.table("popmap.txt",stringsAsFactors=FALSE)
popmap <- popmap[(which(popmap[,1] %in% row.names(snp_geno))),]

pop_designations <- unlist(lapply(1:length(row.names(snp_geno)),function(x) {popmap[which(popmap[,1] %in% row.names(snp_geno)[x]),2]}))

snp_geno <- cbind(pop_designations,snp_geno,stringsAsFactors=FALSE)

library(hierfstat)

if (length(unique(snp_geno[,1]))==1) {
  snp_geno[,1] <- 1
}  
  
snp_geno_stats <- basic.stats(snp_geno)

Ho <- snp_geno_stats$Ho
Hs <- snp_geno_stats$Hs
Fis <- snp_geno_stats$Fis
boot_Fis <- boot.ppfis(snp_geno,nboot=1000,quant=c(0.025,0.975),diploid=TRUE,dig=4)
 
results_matrix <- matrix("",ncol=(dim(Ho)[2]+1),nrow=6)
results_matrix[1,2:(dim(results_matrix)[2])] <- colnames(Ho)
results_matrix[2:(dim(results_matrix)[1]),1] <- c("Ho","He","Fis","lower_Fis_CI","higher_Fis_CI")

for (i in 1:dim(Ho)[2]) {
results_matrix[2,(i+1)] <- mean(Ho[,i],na.rm=TRUE)
results_matrix[3,(i+1)] <- mean(Hs[,i],na.rm=TRUE)
results_matrix[4,(i+1)] <- 1-(mean(Ho[,i],na.rm=TRUE)/mean(Hs[,i],na.rm=TRUE))
}

results_matrix[5:6,(2:dim(results_matrix)[2])] <- t(boot_Fis$fis.ci)
  
basename <- gsub(".vcf","",vcf_name)

print(paste(basename,".He.txt has been written out",sep=""))
write.table(results_matrix,(paste(basename,".He.txt",sep="")),quote=FALSE,row.names=FALSE,col.names=FALSE)  

}  
