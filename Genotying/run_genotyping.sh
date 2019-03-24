#!/bin/bash -l
#SBATCH -A b2010003
#SBATCH -p core -n 8
#SBATCH -J sep vcf
#SBATCH -t 240:00:00
#SBATCH -o /proj/b2014034/private/reseq_analysis/runs/genotyping.output
#SBATCH -e /proj/b2014034/private/reseq_analysis/runs/genotyping.error
#SBATCH --mail-user anna.johansson@scilifelab.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load samtools/1.3
module load bwa/0.7.12
module add FastQC/0.11.2 
module add cutadapt/1.8.0 
module add TrimGalore/0.4.0
module add vcftools

#java -Xmx128g -jar /sw/apps/bioinfo/GATK/3.2.2/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /proj/b2014034/private/reseq_analysis/assembly_mp/N.Backstrom_leptidea.scf.1.4.fasta \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/10.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/11.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/12.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/13.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/14.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/15.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/17.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/18.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/1.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/20.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/21.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/22.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/23.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/24.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/26.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/27.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/29.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/2.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/30.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/31.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/32.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/33.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/34.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/35.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/36.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/37.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/38.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/39.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/3.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/40.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/41.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/42.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/43.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/44.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/46.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/47.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/6.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/7.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/8.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/9.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/Ire-juv-1C.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/Ire-juv-21C.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/Ire-juv-22C.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/Ire-juv-2C.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/Ire-juv-41C.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/Ire-juv-42C.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/Ire-juv-61C.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/Ire-juv-62C.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/Ire-juv-81C.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/Ire-juv-82C.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/Swe-sin-101C.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/Swe-sin-102C.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/Swe-sin-1C.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/Swe-sin-2C.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/Swe-sin-31C.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/Swe-sin-32C.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/Swe-sin-61C.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/Swe-sin-62C.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/Swe-sin-91C.recal.bam.g.vcf \
#-V /proj/b2014034/private/reseq_analysis/bam/gvcf/Swe-sin-92C.recal.bam.g.vcf \
#-o /proj/b2014034/private/reseq_analysis/results/vcf2/genotype_all_20160315.vcf

#java -Xmx128g -jar /sw/apps/bioinfo/GATK/3.2.2/GenomeAnalysisTK.jar \
#-R /proj/b2014034/private/reseq_analysis/assembly_mp/N.Backstrom_leptidea.scf.1.4.fasta \
#-T SelectVariants \
#-V /proj/b2014034/private/reseq_analysis/results/vcf2/genotype_all_20160315.vcf \
#-o /proj/b2014034/private/reseq_analysis/results/vcf2/raw_indels_all_20160315.vcf \
#-selectType INDEL  

#java -Xmx128g -jar /sw/apps/bioinfo/GATK/3.2.2/GenomeAnalysisTK.jar \
#-R /proj/b2014034/private/reseq_analysis/assembly_mp/N.Backstrom_leptidea.scf.1.4.fasta \
#-T SelectVariants -V /proj/b2014034/private/reseq_analysis/results/vcf2/genotype_all_20160315.vcf \
#-o /proj/b2014034/private/reseq_analysis/results/vcf2/raw_snps_all_20160315.vcf \
#-selectType SNP

#minQ_indels=`perl pl/find_qual.pl /proj/b2014034/private/reseq_analysis/results/vcf2/raw_indels_all_20160315.vcf`
#minQ_snps=`perl pl/find_qual.pl /proj/b2014034/private/reseq_analysis/results/vcf2/raw_snps_all_20160315.vcf`

#vcftools --vcf /proj/b2014034/private/reseq_analysis/results/vcf2/raw_indels_all_20160315.vcf --minQ $minQ_indels --recode --out /proj/b2014034/private/reseq_analysis/results/vcf2/raw_indels_all_20160315_qual_filt
#vcftools --vcf /proj/b2014034/private/reseq_analysis/results/vcf2/raw_snps_all_20160315.vcf --minQ $minQ_snps --recode --out /proj/b2014034/private/reseq_analysis/results/vcf2/raw_snps_all_20160315_qual_filt                                                                                                         
                       
#java -Xmx128g -jar /sw/apps/bioinfo/GATK/3.2.2/GenomeAnalysisTK.jar \
#-T VariantRecalibrator \
#-R /proj/b2014034/private/reseq_analysis/assembly_mp/N.Backstrom_leptidea.scf.1.4.fasta \
#-nt 4 \
#-mode SNP \
#-input /proj/b2014034/private/reseq_analysis/results/vcf2/raw_snps_all_20160315.vcf \
#-resource:known=false,training=true,truth=true,prior=15.0 /proj/b2014034/private/reseq_analysis/results/vcf2/raw_snps_all_20160315_qual_filt.recode.vcf \
#-recalFile /proj/b2014034/private/reseq_analysis/bam/recal/snp.tranches.recal \
#-tranchesFile /proj/b2014034/private/reseq_analysis/bam/recal/snp.tranches \
#-an QD \
#-an MQ \
#-an MQRankSum \
#-an ReadPosRankSum \
#-an FS \
#-an DP \
#-tranche 100.0 \
#-tranche 99.9 \
#-tranche 99.0 \
#-tranche 90.0 \
#-rscriptFile /proj/b2014034/private/reseq_analysis/bam/recal/snp.vqsr.r \
#-allPoly  

#java -Xmx128g -jar /sw/apps/bioinfo/GATK/3.2.2/GenomeAnalysisTK.jar \
#-T ApplyRecalibration \
#-R /proj/b2014034/private/reseq_analysis/assembly_mp/N.Backstrom_leptidea.scf.1.4.fasta \
#-input /proj/b2014034/private/reseq_analysis/results/vcf2/raw_snps_all_20160315.vcf \
#-recalFile /proj/b2014034/private/reseq_analysis/bam/recal/snp.tranches.recal \
#-tranchesFile /proj/b2014034/private/reseq_analysis/bam/recal/snp.tranches \
#-o /proj/b2014034/private/reseq_analysis/results/vcf2/recal_snps_all_20160315.vcf \
#-mode SNP  

java -Xmx128g -jar /sw/apps/bioinfo/GATK/3.2.2/GenomeAnalysisTK.jar -T VariantRecalibrator \
-R /proj/b2014034/private/reseq_analysis/assembly_mp/N.Backstrom_leptidea.scf.1.4.fasta \
-nt 4 \
-mode INDEL \
-mG 4 \
-std 10.0 \
-input /proj/b2014034/private/reseq_analysis/results/vcf2/raw_indels_all_20160315.vcf \
-resource:known=false,training=true,truth=true,prior=12.0 /proj/b2014034/private/reseq_analysis/results/vcf2/raw_indels_all_20160315_qual_filt.recode.vcf  \
-recalFile /proj/b2014034/private/reseq_analysis/bam/recal/indel.tranches.recal \
-tranchesFile /proj/b2014034/private/reseq_analysis/bam/recal/indel.tranches \
-an QD \
-an DP \
-an FS \
-an ReadPosRankSum \
-an MQRankSum \
-tranche 100.0 \
-tranche 99.9 \
-tranche 99.9 \
-tranche 90.0 \
-rscriptFile /proj/b2014034/private/reseq_analysis/bam/recal/indel.vqsr.r \
-allPoly  

java -Xmx128g -jar /sw/apps/bioinfo/GATK/3.2.2/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R /proj/b2014034/private/reseq_analysis/assembly_mp/N.Backstrom_leptidea.scf.1.4.fasta \
-input /proj/b2014034/private/reseq_analysis/results/vcf2/raw_indels_all_20160315.vcf \
-recalFile /proj/b2014034/private/reseq_analysis/bam/recal/indel.tranches.recal \
-tranchesFile /proj/b2014034/private/reseq_analysis/bam/recal/indel.tranches \
-o /proj/b2014034/private/reseq_analysis/results/vcf2/recal_indels_all_20160315.vcf \
-mode INDEL

echo Job is done
date

