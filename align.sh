#!/bin/bash

#$ -cwd
#$ -j y
#$ -l m_mem_free=4G
#$ -l tmp_free=400G
#$ -pe threads 8

#### parallel job options.
SAMPLENUMBER LINE
LOGDIR LINE

# Define commonly used files/paths
GATK=$HOME/bin/GenomeAnalysisTK.jar
picard=$HOME/bin/picard-tools-1.126/bin/picard.jar
trimmomatic=$HOME/bin/Trimmomatic-0.32/trimmomatic-0.32.jar

# Get parameters...
PARAMETERFILE LINE
#Parameters=$(sed -n -e "$SGE_TASK_ID p" ParameterFile)
read1=$( echo "$Parameters" | awk '{print $1}' )
read2=$( echo "$Parameters" | awk '{print $2}' )
sample=$( echo "$Parameters" | awk '{print $3}' )
library=$( echo "$Parameters" | awk '{print $4}' )
ref=$( echo "$Parameters" | awk '{print $5}' )
#/sonas-hs/lippman/hpc/home/zlemmon/indexes/SL3.0/SL3.0
#/sonas-hs/lippman/hpc/home/zlemmon/indexes/Pepper155_Genome.fa

######################
# Can just use the samples array to access data later by using standard variable concatenation with strings, but later just use "${samples[@]}". NOTE. If multiple samples (lanes of sequencing) were done on a library separate the samples by semi-colons in the parameters file.
IFS=';' read -r -a read1s <<< "$read1"
IFS=';' read -r -a read2s <<< "$read2"
IFS=';' read -r -a samples <<< "$sample"

#######################
# Rename reads to "$sample"_R[12].fastq.gz
arrlen=$( expr ${#samples[@]} - 1 )
for item in $( seq 0 $arrlen ) ; do
	mv -v "${read1s[$item]}" $TMPDIR/"${samples[$item]}"_R1.fastq.gz
	mv -v "${read2s[$item]}" $TMPDIR/"${samples[$item]}"_R2.fastq.gz
done

##########################
# Run fastqc

echo -ne "\n#######################################################\n"
echo -ne "Beginning fastqc for "$library"\n\n"

mkdir -p ./fastqc
for item in "${samples[@]}" ; do
	/sonas-hs/lippman/hpc/home/zlemmon/bin/FastQC/fastqc -o ./fastqc $TMPDIR/"$item"_R1.fastq.gz
	/sonas-hs/lippman/hpc/home/zlemmon/bin/FastQC/fastqc -o ./fastqc $TMPDIR/"$item"_R2.fastq.gz
done

# start up the analysis of reads...
echo -ne "\n#######################################################\n"
echo -ne "Beginning analysis for "$library"\n" ; date

# align with bwa, and count basic stats with samtools flagstat. Output will be $library_sorted.bam file that has separate alignments for use in a picard MarkDuplicates command and later parallel/chromosome mpileup command.
for item in "${samples[@]}" ; do
	echo -ne "\n###############################################################################################\n"
	echo "Processing sample:$item from library:$library" ; date
	
	##########
	# IMPORTANT FOR PICARD MarkDuplicates
	##########
	# The policy of ENA to have the "Sample_ReadNumber" in front of the actual machine readname causes problems for detecting optical duplicates. uncompressing and awking the readname fixes the problem.
	#echo "Modifying read names..."
	#gunzip -c "$TMPDIR"/"$read1gz" | awk 'NR % 4 == 1 {gsub(/^@ERR[0-9]*.[0-9]* /,"@",$0);print $0;next}{print}' | gzip > "$TMPDIR"/"$read1rn"
	#gunzip -c "$TMPDIR"/"$read2gz" | awk 'NR % 4 == 1 {gsub(/^@ERR[0-9]*.[0-9]* /,"@",$0);print $0;next}{print}' | gzip > "$TMPDIR"/"$read2rn"
	#rm -fv "$TMPDIR"/"$read1gz" "$TMPDIR"/"$read2gz" # remove original reads.

	# Trim reads
	java -Xmx16g -jar "$trimmomatic" PE -threads 8 $TMPDIR/"$item"_R1.fastq.gz $TMPDIR/"$item"_R2.fastq.gz "$TMPDIR"/"$item"_"$library"_P1.fastq.gz "$TMPDIR"/"$item"_"$library"_U1.fastq.gz "$TMPDIR"/"$item"_"$library"_P2.fastq.gz "$TMPDIR"/"$item"_"$library"_U2.fastq.gz ILLUMINACLIP:$HOME/bin/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa:2:40:15:1:FALSE LEADING:30 TRAILING:30 MINLEN:75 TOPHRED33
	# clean out unpaired reads from R1 & R2
	rm -fv "$TMPDIR"/"$item"_"$library"_U1.fastq.gz "$TMPDIR"/"$item"_"$library"_U2.fastq.gz
	# clean out the raw reads from this directory
	rm -fv $TMPDIR/"$item"_R1.fastq.gz $TMPDIR/"$item"_R2.fastq.gz

	# Align library reads (PE) to $ref genome.
	echo "Aligning sample:$item from library:$library (PE) to $ref" ; date
	bwa mem -t 8 -M -R "@RG\tID:$item$library\tSM:$item\tLB:$library\tPL:illumina" "$ref" "$TMPDIR"/"$item"_"$library"_P1.fastq.gz "$TMPDIR"/"$item"_"$library"_P2.fastq.gz > "$TMPDIR"/"$item"_"$library".sam
	rm -fv "$TMPDIR"/"$item"_"$library"_P1.fastq.gz "$TMPDIR"/"$item"_"$library"_P2.fastq.gz
	echo "Initial BWA alignment for sample:$item from library:$library complete" ; date

	# convert to bam, sort, and index
	echo "converting to bam, sorting bam, and indexing bam for sample:$item from library:$library" ; date
	samtools view -b -T $ref -o "$TMPDIR"/"$item"_"$library".bam "$TMPDIR"/"$item"_"$library".sam
	samtools sort -O bam -o "$TMPDIR"/"$item"_"$library"_sorted.bam -T $TMPDIR/temp "$TMPDIR"/"$item"_"$library".bam
	samtools index "$TMPDIR"/"$item"_"$library"_sorted.bam
	# Clean out sam and unsorted bam
	rm -f "$TMPDIR"/"$item"_"$library".sam "$TMPDIR"/"$item"_"$library".bam
	echo "Completed bam conversion and sorting for sample:$item from library:$library" ; date

	echo "Calculating metrics for alignment of sample:$item from library:$library" ; date
	samtools flagstat "$TMPDIR"/"$item"_"$library"_sorted.bam > "$item"_"$library".alignment_metrics
	cat "$item"_"$library".alignment_metrics

	echo "Finished aligning sample:$item from library:$library" ; date
done

echo -ne "\n\n###############################################################################################\n\n"


# For read groups from the same library, mark duplicates with picard and merge bam files.
echo "Marking duplicate reads for library:$library with picard" ; date
# There will be a WARNING thrown for the READ_NAME_REGEX not matching if using default ENA readnames. This is due to the ENA policy to output a "SAMPLENAME.######" Prior to the machine read name "FLOWCELL:LANE:Xcord:Ycord...". The machine read is still there so this problem can be avoided by modifying the readnames prior to alignment and trimming with a simple awk statement.
java -Xmx16g -jar "$picard" MarkDuplicates VALIDATION_STRINGENCY=LENIENT $(printf 'INPUT=%s ' "$TMPDIR"/*"$library"_sorted.bam) OUTPUT="$library"_dups.bam METRICS_FILE="$library"_dups.Metrics VERBOSITY=WARNING 
## index the duplicates marked bam files
samtools index "$library"_dups.bam
echo -ne 'Done!!'

echo -ne "\n\n###############################################################################################\n\n"

echo "Done with primary fastqc, read alignment, and marking duplicates." ; date
