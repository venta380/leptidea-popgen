#!/bin/bash -l
#SBATCH -A b2014034
#SBATCH -p core -n 8
#SBATCH -J vcf_analysis
#SBATCH -t 24:00:00
#SBATCH -o /proj/b2014034/private/reseq_analysis/runs/vcf_analysis.output
#SBATCH -e /proj/b2014034/private/reseq_analysis/runs/vcf_analysis.error
#SBATCH --mail-user anna.johansson@scilifelab.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
module add vcftools
module add samtools
module add BEDTools

######
## Remove indels from vcf files for the three different calling methods, haplotypecaller, freebays and samtools
######

vcftools --gzvcf /proj/b2014034/private/reseq_analysis/results/vcf/genotype_all_20160208.vcf.gz --remove-indels --recode --recode-INFO-all --out /proj/b2014034/private/reseq_analysis/results/vcf/genotype_all_20160208.snps
vcftools --vcf /proj/b2014034/private/reseq_analysis/results/vcf/freebayes.vcf --remove-indels --recode --recode-INFO-all --out /proj/b2014034/private/reseq_analysis/results/vcf/freebayes.snps
vcftools --vcf /proj/b2014034/private/reseq_analysis/results/vcf/samtools_20160222.vcf --remove-indels --recode --recode-INFO-all --out /proj/b2014034/private/reseq_analysis/results/vcf/samtools.snps

######
## Extract the singletons in each fo of the three cases
######

vcftools --vcf /proj/b2014034/private/reseq_analysis/results/vcf/freebayes.snps.recode.vcf --singletons --out /proj/b2014034/private/reseq_analysis/results/vcf/freebayes
vcftools --vcf /proj/b2014034/private/reseq_analysis/results/vcf/samtools.snps.recode.vcf --singletons --out /proj/b2014034/private/reseq_analysis/results/vcf/samtools
vcftools --vcf /proj/b2014034/private/reseq_analysis/results/vcf/genotype_all_20160208.snps.recode.vcf --singletons --out /proj/b2014034/private/reseq_analysis/results/vcf/genotype_all_20160208

######
## Remove singletons from vcf files using perl script
######

perl remove_singletons.pl freebayes freebayes.snps.recode
perl remove_singletons.pl samtools samtools.snps.recode
perl remove_singletons.pl genotype_all_20160208 genotype_all_20160208.snps.recode

bgzip /proj/b2014034/private/reseq_analysis/results/vcf/freebayes.snps.recode.no_singletons.vcf
tabix /proj/b2014034/private/reseq_analysis/results/vcf/freebayes.snps.recode.no_singletons.vcf.gz

bgzip /proj/b2014034/private/reseq_analysis/results/vcf/samtools.snps.recode.no_singletons.vcf 
tabix /proj/b2014034/private/reseq_analysis/results/vcf/samtools.snps.recode.no_singletons.vcf.gz

### Fix vcf-format

cat /proj/b2014034/private/reseq_analysis/results/vcf/genotype_all_20160208.snps.recode.no_singletons.vcf |vcf-convert -v 4.2 > out.vcf
mv out.vcf /proj/b2014034/private/reseq_analysis/results/vcf/genotype_all_20160208.snps.recode.no_singletons.vcf

bgzip /proj/b2014034/private/reseq_analysis/results/vcf/genotype_all_20160208.snps.recode.no_singletons.vcf
tabix /proj/b2014034/private/reseq_analysis/results/vcf/genotype_all_20160208.snps.recode.no_singletons.vcf.gz

#######
## Intersect the methods using vcf-tools. Keep variants present in at least two cases
#######

vcf-isec -f -n +2 /proj/b2014034/private/reseq_analysis/results/vcf/freebayes.snps.recode.no_singletons.vcf.gz /proj/b2014034/private/reseq_analysis/results/vcf/samtools.snps.recode.no_singletons.vcf.gz /proj/b2014034/private/reseq_analysis/results/vcf/genotype_all_20160208.snps.recode.no_singletons.vcf.gz > /proj/b2014034/private/reseq_analysis/results/vcf/golden_set.vcf  

######
## Filter golden_set.vcf to only contain variants that are homozygous for at least one individual.  
######     

perl filter_snp_set.pl golden_set

######
## Sort golden set of variants
######
java -Xmx64g -jar /sw/apps/bioinfo/picard/1.127/milou/picard.jar SortVcf I=/proj/b2014034/private/reseq_analysis/results/vcf/golden_set.filter.vcf O=/proj/b2014034/private/reseq_analysis/results/vcf/golden_set.filter.sorted.vcf SEQUENCE_DICTIONARY=/proj/b2014034/private/reseq_analysis/assembly_mp/N.Backstrom_leptidea.scf.1.4.dict

