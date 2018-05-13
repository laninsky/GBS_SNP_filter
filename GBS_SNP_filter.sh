filename=`cat GBS_SNP_filter.txt | tail -n 1`

if [ -f $filename ]; then
   echo "File $filename exists."
else
   echo "File $filename does not exist."
fi


grep "##" $filename >> header_row.txt
grep -v "##" >> temp
