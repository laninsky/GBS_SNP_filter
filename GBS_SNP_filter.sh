filename=`cat GBS_SNP_filter.txt | tail -n 1`	
basename=`echo $filename | sed 's/.vcf//g'`	

if [ -f $basename.biallelic.vcf ]; then
else
  grep "##" $filename >> header_row.txt	
fi	

Rscript GBS_SNP_filter.R
