sMLH_MLH <- function(vcf_name) {

  if (!require('inbreedR')) install.packages('inbreedR'); library('inbreedR')
  if (!require('vcfR')) install.packages('vcfR'); library('vcfR')
  if (!require('reshape')) install.packages('reshape'); library('reshape')

# Following instructions for importing and reformatting data at https://cran.r-project.org/web/packages/inbreedR/vignettes/inbreedR_step_by_step.html
  print("Reading in file")
  vcf <- read.vcfR(vcf_name, verbose = FALSE )
  print("Extracting genotypes")
  gt <- extract.gt(vcf)
  gt <- as.data.frame(t(gt), stringsAsFactors = FALSE)
  gt[gt == "."] <- NA
  print("Converting formats")
  snp_geno <- do.call(cbind, apply(gt, 2, function(x) colsplit(x, "/", c("a","b"))))
  snp_genotypes <- inbreedR::convert_raw(snp_geno)
  print("If the following prints as TRUE, your genotypes are formatted correctly")                                   
  check_data(snp_genotypes)
  print("Calculating sMLH and MLH")                                 
  sHet <- sMLH(snp_genotypes)
  Het <- MLH(snp_genotypes)

  basename <- gsub(".vcf","",vcf_name)                                 
                                   
  write.table(cbind(colnames(vcf@gt)[-1],as.vector(sHet)),paste(basename,"sMLH.txt",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
  write.table(cbind(colnames(vcf@gt)[-1],as.vector(Het)),paste(basename,"MLH.txt",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)

}                                   
