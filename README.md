# Leptidea Population genetics project
# Generating the VCF file.
Mapping of samples were done using BWA mem. The .bam files obtained from BWA are then processed using the best practices outlined by GATK. The scripts to run the mapping and genotyping are mentioned in the directory titled 'Genotying'. The FASTQ files are submitted to ENA and can be downloaded using the accession number PRJEB21838. The VCF file can be downloaded from the following link https://www.dropbox.com/s/7n1mitj6yg37njk/removed_repeat_content.vcf.gz?dl=0. 

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
The population genetic statistics are calculated in this section.  All scripts for this are placed in the directory ‘Popgen’
## Pre processing 
To run the popgen analysis unbiased to the coverage we will need to extract the total number of positions covered in a window. The per site coverage for each sample can be obtained using the bellow BEDTools command. 
```
for i in $(ls $bamdir*.final.bam)
do
filename=$(basename -- "$i")
filename="${filename%final.bam*}"
genomeCoverageBed -d -ibam $i  > $filename'genomeCoverage.txt'
echo $filename'genomeCoverage.txt' ' ...............Done...............'
done
```
The split_windows_into_10.py takes fasta file and the list of genomeCoverage files to split up files in to 10 parts. This is used in multithreading to make the script run faster. The script get_N_sites.py gives the number of sites covered in each window. The file FINAL_db_N_sites is the final output that contains the number of sites covered in each window.
```
python split_windows_into_10.py  'ls /bam_covrage_files/*.txt | head -10 ' 'Ref.fa'&
python split_windows_into_10.py  'ls /bam_covrage_files/*.txt | head -20 | tail -10' 'Ref.fa'&
python split_windows_into_10.py  'ls /bam_covrage_files/*.txt | head -30 | tail -10' 'Ref.fa'&
python split_windows_into_10.py  'ls /bam_covrage_files/*.txt | head -40 | tail -10' 'Ref.fa'&
python split_windows_into_10.py  'ls /bam_covrage_files/*.txt | head -50 | tail -3' 'Ref.fa'&
wait 

python get_N_sites_V2.0.py -i sinapis           -o  sinapis        	  -p $i &
python get_N_sites_V2.0.py -i reali             -o  reali        	  -p $i &
python get_N_sites_V2.0.py -i sinapis_juvernica -o  sinapis_juvernica -p $i &
#see get_N_sites.sh for scripts run on all samples
wait
```
## Obtaining the allele frequency and nucleotide diversity files
Allele frequency and nucleotide diversity files are used to calculate Dxy, pi, number of fixed, number of shared and number of private polymorphisms in each window across the genome. vcftools was used to obtain allele frequency and nucleotide diversity files using the bellow commands. The command option --max-missing 1.0 is used to output all the SNPs with no missing data
```
vcftools --gzvcf removed_repeat_content.vcf.gz  --keep .lists/juvernica.txt  --max-missing 1.0 --recode --out juvernica
vcftools  --vcf   juvernica.recode.vcf --site-pi --out juvernica.PI.site	
vcftools  --vcf   juvernica.recode.vcf --freq --out  juvernica
```
## Calculating FST
Fst is calculated in windows of 10kb across the genome using the bellow commands. 
```
vcftools --gzvcf removed_repeat_content.vcf.gz  --keep sinapis_juvernica  --max-missing 1.0 --recode --out juvernica_sinapis
vcftools  --vcf   juvernica_sinapis.recode.vcf --weir-fst-pop .lists/juvernica.txt --weir-fst-pop keep .lists/sinapis.txt  --out juvernica_sinapis	
#the mean FST from this output is later scaled to the number of sites covered in the windows. ((MEAN_FST*10000)/N_sites_covered)
```
## Calculating Dxy and fixed positions. 
The function named ‘dxy_window_function’ In house developed python package ‘personal_popgen.py’ was used to calculate DXY. This script also outputs the fixed, shared and private polymorphisms. The script dxy_populations.py uses the function and allele frequency files to calculate DXY. 
## Results
The results from the population genetic analysis has been placed in the directory ‘Popgen/results’. These files contain values for fst,Dxy,Pi, fixed, shared and private. Plots are created using the file Plots.py
# Reference genome:
The *Leptidea sinapis* reference genome can be found in the ENA under the accession number SAMEA104168055. Download it to map of population samples. 


