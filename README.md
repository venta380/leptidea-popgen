# Leptidea Population genetics project
# Generating the VCF file.
Mapping of samples are done using BWA mem. The .bam files obtained from BWA are then processed using the best practices outlined by GATK. The scripts to run the mapping and genotyping are mentioned in the directory titled 'Genotying'. The FASTQ files are submitted to ENA and can be downloaded using the accession number. The VCF file can be downloaded from the following link https://www.dropbox.com/s/7n1mitj6yg37njk/removed_repeat_content.vcf.gz?dl=0. 

# PCA (principal component analysis)
PCA is done using the tool SNPrelate from Bioconductor package in R. The R script PCA_snps.R takes a VCF file as input. This can be manually edited at line number 7 in the script in the variable "vcf.fn". 
# Population structure
The software Admixture was used to analyses the population structure. Cross validation error rate is checked for each value of K before running the Admixture analysis. The bellow are the scripts used to generate this
```
vcftools --gzvcf ../removed_repeat_content.vcf.gz --keep list --positions sinapis.list --recode --out sinapis
vcftools --gzvcf ../removed_repeat_content.vcf.gz --keep list --positions juvernica.list --recode --out juvernica 
vcftools --gzvcf ../removed_repeat_content.vcf.gz --positions all.list --recode --out all

plink --allow-extra-chr --allow-no-sex --noweb --recode12 --vcf sinapis.recode.vcf --make-bed --geno 0 --keep-allele-order --out sinapis
plink --allow-extra-chr --allow-no-sex --noweb --recode12 --vcf juvernica.recode.vcf --make-bed --geno 0 --keep-allele-order --out juvernica
plink --allow-extra-chr --allow-no-sex --noweb --recode12 --vcf all.recode.vcf --make-bed --geno 0 --keep-allele-order --out all

#Cross vaildation error rate estimations
for K in 1 2 3 4 5 6; \
do admixture --cv juvernica.ped $K | tee log${K}.out; done

#Estimation of addmixture proportions for each value of K
for K in 1 2 3 4 5 6; \
do /home/venkat/bin/admixture_linux-1.3.0/admixture  -B2000  all.geno $K; done

```

# Population genetics
The population genetic statistics are calculated in this section all the scripts for this are mentioned in the directory ‘Popgen’
## Pre processing 
The script split.py

# Reference genome:
The *Leptidea sinapis* reference genome can be found in the ENA under the accession number SAMEA104168055. Download it to use the mapping of population samples. 


