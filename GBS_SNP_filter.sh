#GBS_SNP_filter v1.16
if [ `wc -l GBS_SNP_filter.txt | awk '{print $1}'` -eq 8 ]
then

filename=`cat GBS_SNP_filter.txt | head -n 1`	
basename=`echo $filename | sed 's/.vcf//g'`	

if [ -f $basename.biallelic.vcf ]
then
  tempcheck="No"
else
  grep "##" $filename > header_row.txt;
  headerlineno=`wc -l header_row.txt | awk '{print $1}'`;
  headerlineno=$((headerlineno+1))
  tail -n +$headerlineno $filename > temp;
  tempcheck="Yes"
fi	

Rscript GBS_SNP_filter_HWE.R

minrsq=`head -n 5 GBS_SNP_filter.txt | tail -n 1`
for i in `ls *pop.vcf`; do
    ldname=`echo $i | sed 's/.vcf//g'`
    vcftools --vcf $i --plink --out $ldname.plink
    Rscript GBS_SNP_filter_chrom_modifier.R
    plink --allow-extra-chr --r2 --ld-window-r2 $minrsq --file $ldname.plink --out $ldname
    rm -rf *plink*
    rm -rf *nosex
done

Rscript GBS_SNP_filter_rsq.R

if [ $tempcheck == "Yes" ]
then
  rm temp
fi

else
echo "Your GBS_SNP_filter.txt file does not have 8 lines. Are you missing parameters/a blank line on Line 8?"
fi
