# GBS_SNP_filter v1.04

After you get your ddRADseq/GBS dataset through your favourite pipeline, you might want to further filter the SNPs contained in the vcf before doing downstream analysis. This set of scripts filters for bi-allelicness, one SNP/locus, completeness, HWD and LD, and outputs vcf files for each filtering stage.

More detail: this bash/Rscript pipeline first filters for bi-allelic SNPs (and writes out \*.biallelic.vcf), then filters for one SNP/locus (prioritizing the SNP site found in the most individuals. If this is a tie, then the SNP with the highest average coverage. If this is a tie, then randomly selects a SNP. Writes this to \*.oneSNP.vcf). Following this, SNPs are filtered for completeness (according to the parameters you set in GBS_SNP_filter.txt. A vcf file with the format \*.X_Y.vcf will be written out, where X = the proportion of completeness for loci, Y = the proportion of missing data allowed per individual), then for HWE (\*.X_Y.Z_P.HWE.vcf, where Z = the p-value threshold used as a cut-off to suggest a locus is in HWD, and P = the threshold for the number of populations where HWD/LD could occur before that locus was tossed out), and finally for LD (\*.X_Y.Z_P.HWE.Q.ld.vcf, where Q = the Rsq threshold used to chuck out one out of a pair of loci in LD across P populations).

In addition to the filtered vcf files from each step, there will be a number of other files written out: \*.log (contains the number of SNPs/samples across each \*.vcf file, and which samples/loci were binned during any of the steps); \*.X_Y.Z_P.HWE (contains the HWE p-values by row for each locus/by column for each population. This will only include loci that passed the completeness filter); \*.X_Y.Z_P.HWE.\*.pop.vcf (population specific vcf files used for the ld steps. These include only loci that passed the completeness and HWE filters); \*.X_Y.Z_P.HWE.\*.pop.ld (per population files containing pairs of loci in LD as defined by having a RSq >= Q); \*.X_Y.Z_P.HWE.\*.pop.log (log files from PLINK identifying the pairs of loci in LD); and \*.X_Y.Z_P.HWE.Q.rsq (contains loci removed due to LD, and the locus retained that they were in linkage with).

# required inputs
This pipeline requires you have a vcf file output from your favourite pipeline (e.g. ANGSD, ipyrad, stacks) etc., a file called popmap.txt which contains the population code for each individual, and a parameters file called GBS_SNP_filter.txt. These files are described below. It also requires you to have previously installed the R packages dplyr, readr, and stringr, and for you to have R, vcftools, and PLINK in your path. e.g.
```
export PATH="/software/R/bin:/software/bin/vcftools:/software/PLINK/1.09b3.32/plink:$PATH"
```

# input vcf file
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

Depending on how you generated your vcf file, you might have different ideas of what a "locus" is e.g. in the file above, a de novo assembled RAD/GBS dataset, the locus identifier is in the #CHROM column, and the position of the SNP within that locus is given in the POS column. In the following reference-assembled file:
```
##fileformat=VCFv4.2
##fileDate=20180310
##source="Stacks v1.48"
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Allele Depth">
##FORMAT=<ID=GL,Number=.,Type=Float,Description="Genotype Likelihood">
##INFO=<ID=locori,Number=1,Type=Character,Description="Orientation the corresponding Stacks locus aligns in">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  GFO-1   GFO-10  GFO-11  GFO-12  GFO-14  GFO-15  GFO-18  GFO-19  GFO-2
NC_001606.1     14019   223111_6        A       G       .       PASS    NS=17;AF=0.059;locori=p GT:DP:AD        ./.:0:.,.       ./.:0:.,.
NC_031697.1     109273  376732_24       A       G       .       PASS    NS=47;AF=0.074;locori=m GT:DP:AD        ./.:0:.,.       ./.:0:.,.
NC_031697.1     113922  416173_45       G       C       .       PASS    NS=10;AF=0.050;locori=p GT:DP:AD        ./.:0:.,.       ./.:0:.,.
NC_031697.1     114000  416174_43       T       G       .       PASS    NS=10;AF=0.050;locori=m GT:DP:AD        ./.:0:.,.       ./.:0:.,.
NC_031697.1     114021  416174_22       C       G       .       PASS    NS=10;AF=0.050;locori=m GT:DP:AD        ./.:0:.,.       ./.:0:.,.
NC_031697.1     114037  416174_6        G       A       .       PASS    NS=10;AF=0.050;locori=m GT:DP:AD        ./.:0:.,.       ./.:0:.,.
NC_031697.1     121273  468073_5        T       A       .       PASS    NS=2;AF=0.250;locori=p  GT:DP:AD        ./.:0:.,.       ./.:0:.,.
NC_031697.1     123294  408724_12       A       C       .       PASS    NS=11;AF=0.091;locori=m GT:DP:AD        ./.:0:.,.       ./.:0:.,.
NC_031697.1     133307  88175_45        T       C       .       PASS    NS=7;AF=0.071;locori=p  GT:DP:AD        ./.:0:.,.       0/0:3:3,0
NC_031697.1     141382  223355_61       A       T       .       PASS    NS=27;AF=0.093;locori=p GT:DP:AD        ./.:0:.,.       ./.:0:.,.
NC_031697.1     147737  88272_58        C       A       .       PASS    NS=12;AF=0.083;locori=m GT:DP:AD        ./.:0:.,.       0/0:3:3,0
NC_031697.1     147748  88272_47        C       A       .       PASS    NS=8;AF=0.125;locori=m  GT:DP:AD        ./.:0:.,.       0/0:3:3,0
```
The scaffold the locus was assembled against is given in the #CHROM column, while the locus name is given in the ID column separated by the position within that locus by an underscore. For the GBS_SNP_filter.txt file below, you'll need to know what column you want the locus ID based on, and if there is a regular expression needed to remove everything except the locus ID from the contents of that column e.g. for above `_.*`

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
On the first line you need to give the name of your original vcf file that will be processed. On the second line you should give the proportion of samples that has to be equalled or exceeded for a SNP to be retained in the dataset (e.g 0.85 = a SNP needs to be found in 85% or more of total samples to be retained). On the third line you need to give the the proportion of missing loci that will be tolerated for individual samples before they will be removed from the dataset (e.g. 0.9 = a sample can have up to 90% missing data before it is removed from the dataset). On the fourth line you need to give the p-value cut-off for determining whether a locus is out of HWE within a population. On the fifth line you need to give the r^2 cut-off for determining whether SNPs are in LD with each other within a population. On the sixth line you need to give the number of populations a locus has to be out of HWE/in LD across in before that locus is discarded. On the seventh line you need to give the column header for where the locus id you wish to use resides. On the eighth line you need to give any regular expression needed to modify the contents of the column specified on Line 7 to give just the locus name.

e.g. for the first vcf file presented above (make sure to leave a blank line on Line 7 if no regex pattern required)
```
robin.vcf
0.85
0.9
0.05
0.5
3
#CHROM

```

e.g. for the reference-guided vcf file presented above
```
batch_1.vcf
0.65
0.9
0.01
0.5
3
ID
_.*
```

# To run GBS_SNP_filter
Make sure GBS_SNP_filter.sh; GBS_SNP_filter_HWE.R; GBS_SNP_filter_rsq.R; and your inpt vcf file, your popmap.txt, and your GBS_SNP_filter.txt files are located in your working directory. Also make sure you have previously installed the R packages dplyr, readr, and stringr, and have R, vcftools, and PLINK in your path. Then on the terminal:
```
bash GBS_SNP_filter.sh
```
If you are submitting through a slurm system, you might need to preface the bash command with srun within your sbatch script e.g.
```
srun bash GBS_SNP_filter.sh
```

# Troubleshooting
If you find that too many SNPs are being discarded based on the SNP completeness filter (e.g. being found in >= 85% of the samples), it may be that you have had a larger-than-expected number of samples fail. I would suggest changing the second line of GBS_SNP_filter.txt to 0.0 and not filter SNPs based on this metric the first time around. Following filtering of the datasets for samples with high levels of missing data, you could then take the output vcf and run it through another round of filtering, bumping this second line back up to a more stringent value (e.g. 0.85)

# Further utilities
In the utilities folder are scripts for divvying out your "final" vcf into population-specific vcf files, and using reshape, vcfR, hierfstat and inbreedR to calculate individual-level (multi-locus heterozygosity) and population-level (He) estimates of heterozygosity.

# Suggested citation
This code was first published in:
TBD

If you could cite the pub, and the scripts as below, that would be lovely:

Alexander, A. 2018. GBS_SNP_filter v1.x.x. Available from https://github.com/laninsky/GBS_SNP_filter

This pipeline also wouldn't be possible without the programs/packages listed below. Please cite them as well:

R Core Team (2017). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.

Hadley Wickham (2018). stringr: Simple, Consistent Wrappers for Common String Operations. R package version X.X.X. https://CRAN.R-project.org/package=stringr

Hadley Wickham, Romain Francois, Lionel Henry and Kirill Müller (2017). dplyr: A Grammar of Data Manipulation. R package version X.X.X. https://CRAN.R-project.org/package=dplyr

Hadley Wickham, Jim Hester and Romain Francois (2017). readr: Read Rectangular Text Data. R package version X.X.X. https://CRAN.R-project.org/package=readr

Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, Bender D, Maller J, Sklar P, de Bakker PIW, Daly MJ & Sham PC (2007). PLINK: a toolset for whole-genome association and population-based linkage analysis. American Journal of Human Genetics, 81.

Petr Danecek, Adam Auton, Goncalo Abecasis, Cornelis A. Albers, Eric Banks, Mark A. DePristo, Robert Handsaker, Gerton Lunter, Gabor Marth, Stephen T. Sherry, Gilean McVean, Richard Durbin and 1000 Genomes Project Analysis Group. 2011. The Variant Call Format and VCFtools. Bioinformatics. http://dx.doi.org/10.1093/bioinformatics/btr330

Knaus BJ and Grünwald NJ (2017). “VCFR: a package to manipulate and visualize variant call format data in R.” Molecular Ecology Resources, 17(1), pp. 44–53. ISSN 757, http://dx.doi.org/10.1111/1755-0998.12549.

Stoffel, M. A., Esser, M., Kardos, M., Humble, E., Nichols, H., David, P., & Hoffman, J. I. (2016). inbreedR: An R package for the analysis of inbreeding based on genetic markers. Methods in Ecology and Evolution.

H. Wickham. Reshaping data with the reshape package. Journal of Statistical Software, 21(12), 2007.

Goudet, J., 2005. Hierfstat, a package for R to compute and test hierarchical F‐statistics. Molecular Ecology Notes, 5(1), pp.184-186.

# Version history

1.04: Modified GBS_SNP_filter_rsq.R in response to https://github.com/laninsky/GBS_SNP_filter/issues/1: I had been assuming that at least one population would end being dropped due to individual samples failing! 

1.0.3: Fixed a bug introduced in 1.0.2 that was not coping with writing out population specific vcfs for calculating LD if any of the populations were completely missing populations. Added the scripts in the utilities folder.

1.0.2: Fix for bug that meant the first sample in the original vcf appeared in all population-specific vcfs used for calculating HWE and rsq. Modified Rscripts so dplyr, readr and stringr packages were installed if they weren't already when first called.

1.0.1: Fix for bug that arose when all samples in a population were removed due to the missing filter. Minor changes to README.md.

1.0.0: Tested on a single data set so far, and working. If any issues on your data set, please log an issue. Scripts are confirmed to work with these versions of software/libraries: VCFtools/0.1.14-gimkl-2017a-Perl-5.24.1, PLINK/1.09b3.32, R/3.5.0-gimkl-2017a, stringr v 1.3.1, readr v 1.1.1, and dplyr v 0.7.4
