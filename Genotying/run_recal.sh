#!/bin/bash
#

ref=/proj/b2014034/private/reseq_analysis/bam/refs/N.Backstrom_leptidea.scf.1.4.fasta
bwa_ref=/proj/b2014034/private/reseq_analysis/bam/refs/N.Backstrom_leptidea.scf.1.4.fasta
BWA_DIR=/proj/b2014034/private/reseq_analysis/bam/bwa
BAM_DIR=/proj/b2014034/private/reseq_analysis/bam
GVCF_DIR=/proj/b2014034/private/reseq_analysis/bam/gvcf
RUN_DIR=/proj/b2014034/private/reseq_analysis/runs

type=core
cores=8
mem=$(($cores*8))g
account=b2014034
echo $mem

######

for file in $BAM_DIR/*recal.bam;
do

f=$(basename $file);
sample=${f%%_*};       
name=${f%.recal.bam};
out_bam=$BAM_DIR"/"$name".recal.bam";
tmp_bam=$BAM_DIR"/"$name".recal.tmp.bam";
gvcf=$GVCF_DIR"/"$name".recal.bam.g.vcf";
golden_set="/proj/b2014034/private/reseq_analysis/results/vcf/golden_set.filter.sorted.vcf"

run=$RUN_DIR"/recal_"$sample".sh";

echo "#!/bin/bash -l" > $run;
echo "#SBATCH -A $account" >> $run;
echo "#SBATCH -p $type -n $cores" >> $run;
echo "#SBATCH -J recal-$sample" >> $run;
echo "#SBATCH -t 200:00:00" >> $run;
echo "#SBATCH -o /proj/b2014034/private/reseq_analysis/runs/recal_$sample.output" >> $run;
echo "#SBATCH -e /proj/b2014034/private/reseq_analysis/runs/recal_$sample.error" >> $run;
echo "#SBATCH --mail-user anna.johansson@scilifelab.se" >> $run;
echo "#SBATCH --mail-type=ALL" >> $run;

echo "set -euo pipefail" >> $run;

# load some modules

echo "module load bioinfo-tools" >> $run;
echo "module load samtools/1.2" >> $run;
echo "module load bwa/0.7.12" >> $run;

# Step 1 - We need a list of known sites otherwise, GATK will think all the real SNPs in our data are errors.
echo "#java -Xmx$mem -jar /sw/apps/bioinfo/GATK/3.4-46/GenomeAnalysisTK.jar -T BaseRecalibrator -I $file -R $ref -o $BWA_DIR/$name.calibration.csv -knownSites $golden_set" >> $run;
# Step2 - 
echo "#java -Xmx$mem -jar /sw/apps/bioinfo/GATK/3.4-46/GenomeAnalysisTK.jar -T PrintReads -I $file -R $ref -BQSR $BWA_DIR/$name.calibration.csv -o $out_bam" >> $run;

echo "samtools index $out_bam" >> $run;

echo "java -Xmx64g -jar /sw/apps/bioinfo/picard/1.127/milou/picard.jar AddOrReplaceReadGroups I=$out_bam O=$tmp_bam RGID=$name RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$name" >> $run;

echo "mv $tmp_bam $out_bam" >> $run;

echo "samtools index $out_bam" >> $run;

## Variant calling for each sample using HaplotypeCaller, output .g.vcf file
echo "java -Xmx$mem -jar /sw/apps/bioinfo/GATK/3.2.2/GenomeAnalysisTK.jar -T HaplotypeCaller -R $ref -I $out_bam --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o $gvcf" >> $run;

#let me know it is done
#echo "rm $BWA_DIR/*.bam" >> $run;
#echo "rm $BWA_DIR/*.csv" >> $run;
#echo "rm $BWA_DIR/*.intervals" >> $run;
#echo "rm $BWA_DIR/*.metrics" >> $run;

echo "echo "Job is done"" >> $run;
echo "date" >> $run;

#sbatch $run;

done
