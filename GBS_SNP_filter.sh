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
for i in `ls *pop.vcf`; do
    (
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
