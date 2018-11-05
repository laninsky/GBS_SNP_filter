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
The sMLH results will be suffixed by "sMLH.txt", and MLH results suffixed by "MLH.txt"

# pop_He
Given a vcf file, this script will calculate the gene diversity, or heterozygosity (He), of the population for each locus, and the effective number of alleles in the population (through vcfR). You will also need to specify the minimum number of individuals you require at a locus for it to be output. You call the script (after pasting the whole thing into R and/or sourcing it) by:
```
pop_He(vcf_name,numberofheaders)
# e.g.
pop_He("hihi.biallelic.Bushy_Park.pop.vcf",5)
```
