#!/bin/bash -l
#SBATCH -A b2014034
#SBATCH -p core -n 1 
#SBATCH -J samtools
#SBATCH -t 180:00:00
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

samtools mpileup -ugf /proj/b2014034/private/reseq_analysis/assembly_mp/N.Backstrom_leptidea.scf.1.4.fasta \
/proj/b2014034/private/reseq_analysis/bam/10.final.bam \
/proj/b2014034/private/reseq_analysis/bam/11.final.bam \
/proj/b2014034/private/reseq_analysis/bam/12.final.bam \
/proj/b2014034/private/reseq_analysis/bam/13.final.bam \
/proj/b2014034/private/reseq_analysis/bam/14.final.bam \
/proj/b2014034/private/reseq_analysis/bam/15.final.bam \
/proj/b2014034/private/reseq_analysis/bam/17.final.bam \
/proj/b2014034/private/reseq_analysis/bam/18.final.bam \
/proj/b2014034/private/reseq_analysis/bam/1.final.bam \
/proj/b2014034/private/reseq_analysis/bam/20.final.bam \
/proj/b2014034/private/reseq_analysis/bam/21.final.bam \
/proj/b2014034/private/reseq_analysis/bam/22.final.bam \
/proj/b2014034/private/reseq_analysis/bam/23.final.bam \
/proj/b2014034/private/reseq_analysis/bam/24.final.bam \
/proj/b2014034/private/reseq_analysis/bam/26.final.bam \
/proj/b2014034/private/reseq_analysis/bam/27.final.bam \
/proj/b2014034/private/reseq_analysis/bam/29.final.bam \
/proj/b2014034/private/reseq_analysis/bam/30.final.bam \
/proj/b2014034/private/reseq_analysis/bam/31.final.bam \
/proj/b2014034/private/reseq_analysis/bam/32.final.bam \
/proj/b2014034/private/reseq_analysis/bam/33.final.bam \
/proj/b2014034/private/reseq_analysis/bam/34.final.bam \
/proj/b2014034/private/reseq_analysis/bam/35.final.bam \
/proj/b2014034/private/reseq_analysis/bam/36.final.bam \
/proj/b2014034/private/reseq_analysis/bam/37.final.bam \
/proj/b2014034/private/reseq_analysis/bam/38.final.bam \
/proj/b2014034/private/reseq_analysis/bam/39.final.bam \
/proj/b2014034/private/reseq_analysis/bam/3.final.bam \
/proj/b2014034/private/reseq_analysis/bam/40.final.bam \
/proj/b2014034/private/reseq_analysis/bam/41.final.bam \
/proj/b2014034/private/reseq_analysis/bam/42.final.bam \
/proj/b2014034/private/reseq_analysis/bam/43.final.bam \
/proj/b2014034/private/reseq_analysis/bam/44.final.bam \
/proj/b2014034/private/reseq_analysis/bam/46.final.bam \
/proj/b2014034/private/reseq_analysis/bam/47.final.bam \
/proj/b2014034/private/reseq_analysis/bam/6.final.bam \
/proj/b2014034/private/reseq_analysis/bam/7.final.bam \
/proj/b2014034/private/reseq_analysis/bam/8.final.bam \
/proj/b2014034/private/reseq_analysis/bam/9.final.bam \
/proj/b2014034/private/reseq_analysis/bam/Ire-juv-1C.final.bam \
/proj/b2014034/private/reseq_analysis/bam/Ire-juv-21C.final.bam \
/proj/b2014034/private/reseq_analysis/bam/Ire-juv-22C.final.bam \
/proj/b2014034/private/reseq_analysis/bam/Ire-juv-2C.final.bam \
/proj/b2014034/private/reseq_analysis/bam/Ire-juv-41C.final.bam \
/proj/b2014034/private/reseq_analysis/bam/Ire-juv-42C.final.bam \
/proj/b2014034/private/reseq_analysis/bam/Ire-juv-61C.final.bam \
/proj/b2014034/private/reseq_analysis/bam/Ire-juv-62C.final.bam \
/proj/b2014034/private/reseq_analysis/bam/Ire-juv-81C.final.bam \
/proj/b2014034/private/reseq_analysis/bam/Ire-juv-82C.final.bam \
/proj/b2014034/private/reseq_analysis/bam/Swe-sin-101C.final.bam \
/proj/b2014034/private/reseq_analysis/bam/Swe-sin-102C.final.bam \
/proj/b2014034/private/reseq_analysis/bam/Swe-sin-1C.final.bam \
/proj/b2014034/private/reseq_analysis/bam/Swe-sin-2C.final.bam \
/proj/b2014034/private/reseq_analysis/bam/Swe-sin-31C.final.bam \
/proj/b2014034/private/reseq_analysis/bam/Swe-sin-32C.final.bam \
/proj/b2014034/private/reseq_analysis/bam/Swe-sin-61C.final.bam \
/proj/b2014034/private/reseq_analysis/bam/Swe-sin-62C.final.bam \
/proj/b2014034/private/reseq_analysis/bam/Swe-sin-91C.final.bam \
/proj/b2014034/private/reseq_analysis/bam/Swe-sin-92C.final.bam | bcftools call -vmO z -o /proj/b2014034/private/reseq_analysis/results/vcf/samtools_20160222.vcf.gz

echo Job is done
date
