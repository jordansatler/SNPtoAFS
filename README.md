# SNPtoAFS

This code will convert a SNP data set into a folded joint allele frequency
spectrum (AFS) for use in *fastsimcoal2* (Excoffier et al. 2013). This will
convert SNP data to an AFS, using either linked or unlinked SNPs, and will
allow for missing data with the user specifying a threshold for the amount
of missing data in both populations (for example, 75%). The script will 
remove all SNPs that do not meet the minimum threshold; for those SNPs that 
do, if they meet the threshold, the minor allele will be counted. For those 
SNPs that exceed the threshold in a population, the script with subsample 
with replacement until enough alleles are sampled to meet the threshold, 
then the minor allele will be counted. The use may also specify the number 
of replicates to perform, to account for variation in the downsampling 
procedure. 

Two scripts are included in the repository: 

___
**AFS\_FSC\_total.py** will convert a SNP data set to a folded two-population 
AFS using the minor allele counts for each SNP. This will handle files 
from AftrRAD (Sovic et. al 2015) and pyRAD (Eaton 2014). Users specify
the sampling threshold (as a percentage; 75 is used if you want a 75%
threshold) in both populations, number of replicated AFS, and if you want
the AFS to built with linked SNPs or to subsample a single SNP at random.

**SNPtoAFSready.py** will convert a *.snps file from pyRAD for use with the 
AFS\_FSC\_total.py script. 
___

The traits file gives individuals linked to populations.

Ex:
Individuals	Populations
BL1a	west
BL1b	west
BL2a	west
BL2b	west
AS6a	east
AS6b	east


References:



