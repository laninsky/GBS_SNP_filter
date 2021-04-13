sMLH_MLH <- function(vcf_name,parallel=TRUE,redo=FALSE) {
  
  check<-check_continue(vcf_name,redo=redo)
  if(check==TRUE){
    basename <- gsub(".vcf","",vcf_name)                                 
    
    if (!require('inbreedR')) install.packages('inbreedR'); library('inbreedR')
    if (!require('vcfR')) install.packages('vcfR'); library('vcfR')
    if (!require('reshape')) install.packages('reshape'); library('reshape')
    if (!require('future.apply')) install.packages('future.apply'); library('future.apply')
    
    if(parallel==FALSE){plan(sequential)}else{plan(multiprocess)}
    
    if(length(list.files(pattern=paste(basename,"snp_genotypes.RData",sep="_")))==1){
      print(paste("Reading in previously generated snp genotypes for",vcf_name))
      load(paste(basename,"snp_genotypes.RData",sep="_"))
      print("Also reading vcf file to get sample names")
      vcf <- read.vcfR(vcf_name, verbose = FALSE ) # needed for colnames
    }else{
      # Following instructions for importing and reformatting data at https://cran.r-project.org/web/packages/inbreedR/vignettes/inbreedR_step_by_step.html
      print("Reading in VCF file")
      vcf <- read.vcfR(vcf_name, verbose = FALSE )
      print("Extracting genotypes")
      gt <- extract.gt(vcf)
      gt <- as.data.frame(t(gt), stringsAsFactors = FALSE)
      gt[gt == "."] <- NA
      print("Finding and removing any loci with no data for either allele. This step may take a while.")
      all_NAs <- unlist(future_lapply(1:(dim(gt)[2]),function(x){all(is.na(gt[1:(dim(gt)[1]),x]))})) 
      gt <- gt[,which(all_NAs[]==FALSE)]  
      print("Converting formats. This step will also take a while.")  
      snp_geno <- do.call(cbind, future_apply(gt, 2, function(x) colsplit(x, "/", c("a","b"))))
      snp_genotypes <- inbreedR::convert_raw(snp_geno)
      print("If the following prints as TRUE, your genotypes are formatted correctly")                                   
      print(check_data(snp_genotypes))
      save(snp_genotypes,file=paste(basename,"snp_genotypes.RData",sep="_"))
    }
    
    print("Calculating sMLH and MLH")                                 
    sHet <- sMLH(snp_genotypes)
    Het <- MLH(snp_genotypes)
    
    write.table(cbind(colnames(vcf@gt)[-1],as.vector(sHet)),paste(basename,".sMLH.txt",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
    write.table(cbind(colnames(vcf@gt)[-1],as.vector(Het)),paste(basename,".MLH.txt",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
    
  }
}

check_continue<-function(vcf_name,redo=FALSE){
  basename <- gsub(".vcf","",vcf_name)                                 
  if(length(list.files(pattern=paste(basename,".sMLH.txt",sep="")))==1){
    if(redo==FALSE){
      print("Output files for this vcf already exist. Set redo=TRUE to overwrite them")
      return(FALSE)
    }else{
      print("overwriting")
      return(TRUE)
    }
  }else{
    return(TRUE)
  }
}
