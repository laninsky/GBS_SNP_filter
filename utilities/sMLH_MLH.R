if (!require('inbreedR')) install.packages('inbreedR'); library('inbreedR')
if (!require('vcfR')) install.packages('vcfR'); library('vcfR')
if (!require('reshape')) install.packages('reshape'); library('reshape')

nstall.packages("data.table")
install.packages("inbreedR")
library(inbreedR)

install.packages("ape")
install.packages("vegan")
install.packages("vcfR")
install.packages("reshape")
library(vcfR)
library(reshape)

# Following instructions for importing and reformatting data at https://cran.r-project.org/web/packages/inbreedR/vignettes/inbreedR_step_by_step.html
vcf_file <- "hihi_prev_filter.0.85_0.9.0.05_3.HWE.0.5.ld.vcf"
vcf <- read.vcfR(vcf_file, verbose = FALSE )
gt <- extract.gt(vcf)
gt <- as.data.frame(t(gt), stringsAsFactors = FALSE)
gt[gt == "."] <- NA
snp_geno <- do.call(cbind, apply(gt, 2, function(x) colsplit(x, "/", c("a","b"))))
hihi_snp_genotypes <- inbreedR::convert_raw(snp_geno)
check_data(hihi_snp_genotypes)

sHet <- sMLH(hihi_snp_genotypes)
Het <- MLH(hihi_snp_genotypes)

write.table(cbind(colnames(vcf@gt)[-1],as.vector(sHet)),"sHet_hihi.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

write.table(cbind(colnames(vcf@gt)[-1],as.vector(Het)),"Het_hihi.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
