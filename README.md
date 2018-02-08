# Leptidea Popgen project

# PCA (principal component analysis)
The R script PCA_snps.R takes input as a VCF file. This can be manually edited at line number 7 in the script in the variable "vcf.fn". 
# Population structure
The software Admixture was used to analyse the population structure. The bellow are the scripts used to generate this
```
vcftools --gzvcf ../removed_repeat_content.vcf.gz --keep list --positions sinapis.list --recode --out sinapis
vcftools --gzvcf ../removed_repeat_content.vcf.gz --keep list --positions juvernica.list --recode --out juvernica 
vcftools --gzvcf ../removed_repeat_content.vcf.gz --positions all.list --recode --out all

plink --allow-extra-chr --allow-no-sex --noweb --recode12 --vcf sinapis.recode.vcf --make-bed --geno 0 --keep-allele-order --out sinapis
plink --allow-extra-chr --allow-no-sex --noweb --recode12 --vcf juvernica.recode.vcf --make-bed --geno 0 --keep-allele-order --out juvernica
plink --allow-extra-chr --allow-no-sex --noweb --recode12 --vcf all.recode.vcf --make-bed --geno 0 --keep-allele-order --out all

for K in 2; \
do admixture --cv juvernica.ped $K | tee log${K}.out; done

for K in 1 2 3 4 5 6; \
do /home/venkat/bin/admixture_linux-1.3.0/admixture  -B2000  all.geno $K; done

```
# Reference genome:
# Reference genome:
The leptidea refrence genome can be found in the ENA under the accesion number SAMEA104168055. Download it to use the mapping of population samples. 

