#!/bin/bash
#

####
## Input a vcf file with the 2*10 individuals in the two populations
## ./snp_selection.sh /proj/b2014034/private/reseq_analysis/results/vcf2/recal_snps_all_20160315.vcf_g4_Swe.vcf
###

module add bioinfo-tools vcftools samtools BEDOPS vep/82 BEDTools/2.23.0

ref=/proj/b2014034/private/reseq_analysis/bam/refs/N.Backstrom_leptidea.scf.1.4.fasta
REP=/proj/b2014034/private/reseq_analysis/assembly_mp/repetetive_positions.bed 
GAP=/proj/b2014034/private/reseq_analysis/assembly_mp/gap_positions.bed
vcf=$1

######

if [ ! -f $vcf.rep.vcf ]; then
intersectBed -v -header -a $vcf -b $REP > $vcf.rep.vcf 
fi

if [ ! -f $vcf.rep.vcf.maf.bed ]; then
perl pl/filter_maf.pl $vcf.rep.vcf
fi

if [ ! -f $vcf.rep.vcf.maf.bed.no_gap.bed ]; then
intersectBed -v -a $vcf.rep.vcf.maf.bed -b $GAP > $vcf.rep.vcf.maf.bed.no_gap.bed
fi

if [ ! -f $vcf.list.txt ]; then
perl pl/filter_list.pl $vcf
fi


