# GBS_SNP_filter

First filter for bi-allelic SNPs
Then take one SNP per locus







#Then
filtering SNPs based on LD and HWE status

p1 = SNP 1 allele 1 (A1) freq
p2 = SNP 1 allele 2 (A2) freq
q1 = SNP 2 allele 1 (B1) freq
q2 = SNP 2 allele 2 (B2) freq

observed vs expected
A1B1 freq x11, expected p1*q1
A1B2 freq x12, expected p1*q2
A2B1 freq x21, expected p2*q1
A2B2 freq x22, expected p2*q2

D = (x11)(x22) – (x12)(x21)
D = x11 – p1q1

D' = D/(min(p1q1, p2q2) when D < 0 or 
D' = D/(min(p1q2, p1q2) when D > 0 

r^2 = D^2/p1q1p2q2

current plan: for each population (write out on the fly to variables that have something like "LD_log" or something we can read back in), stick D' below the diagonal and r^2 above the diagonal

output D' by r^2 for original data for each population plot (might help with people making cut-offs)

Do this after filtering for SNPs too.

setting number of pops has to be in LD within, r^2 and D' cut-offs for filtering.

PLINK might be able to do this rather than coding from scratch
https://bioinformatics.stackexchange.com/questions/3040/ld-analysis-in-plink-based-on-reference-and-a-snp-list


https://www.ndsu.edu/pubweb/~mcclean/plsc731/Linkage%20Disequilibrium%20-%20Association%20Mapping%20in%20Plants-lecture-overheads.pdf
