pop_He <- function(vcf_name) {

  if (!require('vcfR')) install.packages('vcfR'); library('vcfR')

vcf_file <- "hihi_prev_filter.0.85_0.9.0.05_3.HWE.0.5.ld.vcf"
vcf <- read.vcfR(vcf_file, verbose = FALSE )
popsum_total.txt <- gt.to.popsum(vcf)
