# GBS_SNP_filter

This bash/Rscript pipeline first filters for bi-allelic SNPs (and writes out \*.biallelic.vcf), then filters for one SNP/locus (prioritizing the SNP found in the most individuals. If this is a tie, then the SNP with the highest average coverage. If this is a tie, then randomly selects a SNP. Writes this to \*.oneSNP.vcf). Following this, SNPs are filtered for completeness (according to the parameters you set in GBS_SNP_filter.txt. A vcf file with the format \*.X_Y.vcf will be written out, where X = the proportion of completeness for loci, Y = the proportion of missing data allowed per individual), then for HWE (\*.X_Y.HWE.vcf), and finally for LD (\*.X_Y.final.vcf)

This pipeline requires you have a vcf file output from your favourite pipeline (e.g. ANGSD, ipyrad, stacks) etc., a file called popmap.txt which contains the population code for each individual, and a parameters file called GBS_SNP_filter.txt. These files are described below. It also requires you to have previously installed the R package tidyverse.

# vcf file
Header lines starting with "##" will be ignored. The script expects your first sample to be in column 10. Samplenames should not have the search term "\_cov" in them as this will be used for filtering coverage within the script.
```
##fileformat=VCFv4.0
##fileDate=2018/05/09
##source=ipyrad_v.0.7.24
##reference=pseudo-reference (most common base at site)
##phasing=unphased
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=CATG,Number=1,Type=String,Description="Base Counts (CATG)">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  B104599-GAAGCATCA       B104930-TGGCAGCA        B104957-GAACCGTAA       B109462-TTGTCAGAA       B109466-GAACCTGA
locus_2 3       .       G       A       13      PASS    NS=21;DP=179    GT:DP:CATG      1/1:9:0,9,0,0   0/0:6:0,0,0,6   ./.:0:0,0,0,0   ./.:0:0,0,0,0   ./.:0:0,0,0,0   0/1:12:0,9,0,3  ./.:0:
locus_20        22      .       C       T       13      PASS    NS=99;DP=1182   GT:DP:CATG      0/0:8:8,0,0,0   0/0:14:14,0,0,0 0/0:9:9,0,0,0   0/0:8:8,0,0,0   0/0:8:8,0,0,0   0/0:15:15,0,0,
locus_21        11      .       C       T       13      PASS    NS=8;DP=65      GT:DP:CATG      0/0:8:8,0,0,0   ./.:0:0,0,0,0   ./.:0:0,0,0,0   ./.:0:0,0,0,0   ./.:0:0,0,0,0   ./.:0:0,0,0,0
locus_21        23      .       C       T       13      PASS    NS=8;DP=65      GT:DP:CATG      0/0:8:8,0,0,0   ./.:0:0,0,0,0   ./.:0:0,0,0,0   ./.:0:0,0,0,0   ./.:0:0,0,0,0   ./.:0:0,0,0,0
locus_21        40      .       G       A       13      PASS    NS=8;DP=65      GT:DP:CATG      0/0:8:0,0,0,8   ./.:0:0,0,0,0   ./.:0:0,0,0,0   ./.:0:0,0,0,0   ./.:0:0,0,0,0   ./.:0:0,0,0,0
locus_29        18      .       A       G       13      PASS    NS=113;DP=1538  GT:DP:CATG      1/1:8:0,0,0,8   0/0:14:0,14,0,0 0/0:16:0,16,0,0 1/0:12:0,6,0,6  0/0:10:0,10,0,0 0/0:16:0,16,0,
locus_29        44      .       T       C       13      PASS    NS=113;DP=1538  GT:DP:CATG      1/1:8:8,0,0,0   0/1:14:10,0,4,0 0/0:16:0,0,16,0 0/1:12:6,0,6,0  1/1:10:10,0,0,0 0/1:16:7,0,9,0
locus_32        28      .       G       A       13      PASS    NS=53;DP=569    GT:DP:CATG      1/1:8:0,8,0,0   ./.:0:0,0,0,0   ./.:0:0,0,0,0   ./.:0:0,0,0,0   0/0:12:0,0,0,12 1/1:10:0,10,0,
locus_33        38      .       C       T       13      PASS    NS=110;DP=1543  GT:DP:CATG      1/1:8:0,0,8,0   1/0:17:6,0,11,0 1/1:10:0,0,10,0 1/1:18:0,0,18,0 1/0:13:5,0,8,0  0/0:16:16,0,0,
locus_35        72      .       C       A       13      PASS    NS=53;DP=72     GT:DP:CATG      0/0:1:1,0,0,0   ./.:0:0,0,0,0   ./.:0:0,0,0,0   ./.:0:0,0,0,0   ./.:0:0,0,0,0   ./.:0:0,0,0,0
locus_35        74      .       G       A       13      PASS    NS=47;DP=56     GT:DP:CATG      0/0:1:0,0,0,1   ./.:0:0,0,0,0   ./.:0:0,0,0,0   ./.:0:0,0,0,0   ./.:0:0,0,0,0   ./.:0:0,0,0,0

```

# popmap.txt
Sample name (exactly matching that in the vcf file) in the left hand column, separated by white space from population name in the right column
```
B104599-GAAGCATCA Kai
B104930-TGGCAGCA  Kai
B104957-GAACCGTAA Kai
B109462-TTGTCAGAA Kai
B109466-GAACCTGA  Nuku
B109467-AAGAGTCG  Nuku
B109469-AACATCGCA Nuku
B109470-CGTGGTGCA Nuku
```

# GBS_SNP_filter.txt
On the first line you need to give the name of your original vcf file that will be processed. On the second line you should give the proportion of samples that has to be equalled or exceeded for a SNP to be retained in the dataset (e.g 0.85 = a SNP needs to be found in 85% or more of total samples to be retained). On the third line you need to give the the proportion of missing loci that will be tolerated for individual samples before they will be removed from the dataset (e.g. 0.9 = a sample can have up to 90% missing data before it is removed from the dataset). On the fourth line you need to give the p-value alpha cut-off for determining whether a locus is out of HWE within a population. On the fifth line you need to give the r^2 cut-off for determining whether SNPs are in LD with each other within a population. On the sixth line you need to give the number of populations a locus has to be out of HWE/in LD across in before that locus is discarded. 
```
robin.vcf
0.85
0.9
0.05
0.5
3
```
