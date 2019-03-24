#!/bin/bash
#

ref=/proj/b2014034/private/reseq_analysis/bam/refs/N.Backstrom_leptidea.scf.1.4.fasta
bwa_ref=/proj/b2014034/private/reseq_analysis/bam/refs/N.Backstrom_leptidea.scf.1.4.fasta
FASTQ_DIR=/proj/b2014034/private/reseq_analysis/trimmed_fastq3
BWA_DIR=/proj/b2014034/private/reseq_analysis/bam/bwa
BAM_DIR=/proj/b2014034/private/reseq_analysis/bam
RUN_DIR=/proj/b2014034/private/reseq_analysis/runs
cores=8
mem=$(($cores*8))g

echo $mem

######

#for file in $FASTQ_DIR/*_R1_001_val_1.fq.gz;
for file in $FASTQ_DIR/Swe-*_R1_001_val_1.fq.gz;
#for file in $FASTQ_DIR/Ire-*_R1_001_val_1.fq.gz;
do

f=$(basename $file);
sample=${f%%_*};       
name=${f%_R1_001_val_1.fq.gz};

run=$RUN_DIR"/mapping_"$sample".sh";

echo "#!/bin/bash -l" > $run;
echo "#SBATCH -A b2010003" >> $run;
echo "#SBATCH -p core -n $cores" >> $run;
echo "#SBATCH -J mapping-$sample" >> $run;
echo "#SBATCH -t 200:00:00" >> $run;
echo "#SBATCH -o /proj/b2014034/private/reseq_analysis/runs/mapping_$sample.output" >> $run;
echo "#SBATCH -e /proj/b2014034/private/reseq_analysis/runs/mapping_$sample.error" >> $run;
echo "#SBATCH --mail-user anna.johansson@scilifelab.se" >> $run;
echo "#SBATCH --mail-type=ALL" >> $run;

# load some modules
echo "module load bioinfo-tools" >> $run;
echo "module load samtools/1.2" >> $run;
echo "module load bwa/0.7.12" >> $run;

fasta1=$FASTQ_DIR/$name"_R1_001_val_1.fq.gz";
fasta2=$FASTQ_DIR/$name"_R2_001_val_2.fq.gz";

# Go from fastQ to BAM #

echo "bwa mem -t $cores -M $bwa_ref -R '@RG\tID:$sample\tSM:$sample\tPL:ILLUMINA' $fasta1 $fasta2 | samtools import $ref.fai - - | samtools sort - $BWA_DIR/$sample" >> $run;
echo "samtools index $BWA_DIR/$sample.bam" >> $run;

# Delete and mark duplicates using PICARD
echo "java -Xmx$mem -jar /sw/apps/bioinfo/picard/1.127/milou/picard.jar MarkDuplicates INPUT=$BWA_DIR/$sample.bam OUTPUT=$BWA_DIR/$sample.bam.dedup.bam METRICS_FILE=$BWA_DIR/$sample.bam.metrics REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT" >> $run

echo "samtools index $BWA_DIR/$sample.bam.dedup.bam" >> $run;

# Process reads with GATK
# Step 1 - realign locally around potential indels. First, identify possible sites to realign.
echo "java -Xmx$mem -jar /sw/apps/bioinfo/GATK/3.4-46/GenomeAnalysisTK.jar -I $BWA_DIR/$sample.bam.dedup.bam -R $ref -T RealignerTargetCreator -o $BWA_DIR/$sample.intervals" >> $run;

# Feed the intervals file back into GATK with a different argument to actually do the realignments
echo "java -Xmx$mem -jar /sw/apps/bioinfo/GATK/3.4-46/GenomeAnalysisTK.jar -I $BWA_DIR/$sample.bam.dedup.bam -R $ref -T IndelRealigner -o $BWA_DIR/$sample.bam.dedup.realign.bam -targetIntervals $BWA_DIR/$sample.intervals" >> $run;

# Reindexing the picard file
echo "samtools index $BWA_DIR/$sample.bam.dedup.realign.bam" >> $run;

# Perform quality recalibration with GATK. #

# Step 1 - We need a list of known sites otherwise, GATK will think all the real SNPs in our data are errors.
echo "java -Xmx$mem -jar /sw/apps/bioinfo/GATK/3.4-46/GenomeAnalysisTK.jar -T BaseRecalibrator -I $BWA_DIR/$sample.bam.dedup.realign.bam -R $ref -o $BWA_DIR/$sample.bam.dedup.realign.calibration.csv -knownSites /proj/b2014034/private/reseq_analysis/assembly_mp/knownSites.vcf" >> $run;
# Step2 - 
echo "java -Xmx$mem -jar /sw/apps/bioinfo/GATK/3.4-46/GenomeAnalysisTK.jar -T PrintReads -I $BWA_DIR/$sample.bam.dedup.realign.bam -R $ref -BQSR $BWA_DIR/$sample.bam.dedup.realign.calibration.csv -o $BAM_DIR/$sample.final.bam" >> $run;
echo "samtools index $BAM_DIR/$sample.final.bam" >> $run;

# Quality test #
# Doing flagstat to get stats

echo "samtools flagstat $BAM_DIR/$sample.final.bam > $BAM_DIR/flagstat/flagstat.$sample.final.bam" >> $run;

echo "qualimap bamqc -bam $BAM_DIR/$sample.final.bam -outdir $BAM_DIR/qualimap/$sample" >> $run;

## Variant calling for each sample using HaplotypeCaller, output .g.vcf file
echo "java -Xmx$mem -jar /sw/apps/bioinfo/GATK/3.2.2/GenomeAnalysisTK.jar -T HaplotypeCaller -R $ref -I $BAM_DIR/$sample.final.bam --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o $BAM_DIR/$sample.final.bam.g.vcf" >> $run

#let me know it is done
#echo "rm $BWA_DIR/*.bam" >> $run;
#echo "rm $BWA_DIR/*.csv" >> $run;
#echo "rm $BWA_DIR/*.intervals" >> $run;
#echo "rm $BWA_DIR/*.metrics" >> $run;

echo "echo "Job is done"" >> $run;
echo "date" >> $run;

#sbatch $run;

done
