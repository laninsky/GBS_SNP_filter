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

Rscript GBS_SNP_filter_HWE.R

minrsq=`head -n 5 GBS_SNP_filter.txt | tail -n 1`
for i in `ls *pop.vcf`; do
    ldname=`echo $i | sed 's/.vcf//g'`
    vcftools --vcf $i --plink --out $ldname.plink
    plink --r2 --ld-window-r2 $minrsq --file $ldname.plink --out $ldname
    rm -rf *plink*
    rm -rf *nosex
done

