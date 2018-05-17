filename=`cat GBS_SNP_filter.txt | head -n 1`	
basename=`echo $filename | sed 's/.vcf//g'`	

if [ -f $basename.biallelic.vcf ]
then
  :
else 
  grep "##" $filename >> header_row.txt;
  headerlineno=`wc -l header_row.txt | awk '{print $1}'`;
  headerlineno=$((headerlineno+1))
  tail -n +$headerlineno $filename > temp;
fi	

Rscript GBS_SNP_filter.R
