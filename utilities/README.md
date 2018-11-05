# split_vcf_by_pop.R
Given a popmap.txt file (see main folder) in the same directory as the vcf file you want to split into population-specific vcf files, this R-script will divvy up your original vcf file into population-specific vcf files with the suffix: "popname.pop.vcf". You'll need to know the number of header lines in your original vcf file (numberofheaders). You call the script (after pasting the whole thing into R and/or sourcing it) by:
```
split_vcf_by_pop(vcf_name,numberofheaders)
# e.g.
split_vcf_by_pop("hihi.biallelic.vcf",10)
```

# sMLH_MLH.R
Given a vcf file, this script will calculate MLH and sMLH for all individuals in your file. You call the script (after pasting the whole thing into R and/or sourcing it) by:
```
sMLH_MLH(vcf_name)
# e.g.
sMLH_MLH("hihi.biallelic.vcf")
```
