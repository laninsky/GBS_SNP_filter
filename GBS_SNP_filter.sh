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

N=`tail -n 1 GBS_SNP_filter.txt`
minrsq=`head -n 5 GBS_SNP_filter.txt | tail -n 1`
for i in `ls *pop.vcf`; do
    (
    ldname=`echo $i | sed 's/.vcf//g'`
    vcftools --vcf $i --geno-r2 --min-r2 $minrsq --out $ldname
  # .. do your stuff here
        echo "starting task $i.."
        sleep $(( (RANDOM % 3) + 1))
    ) &

    # allow only to execute $N jobs in parallel
    if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
        # wait only for first job
        wait -n
    fi

done
wait
