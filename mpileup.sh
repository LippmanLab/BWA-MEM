#!/bin/bash
#$ -cwd 
#$ -j y
#$ -l m_mem_free=2G
#$ -pe threads 14

#$ -t 1
LOGDIR LINE

# Define commonly used files/paths
GATK=$HOME/bin/GenomeAnalysisTK.jar
picard=$HOME/bin/picard-tools-1.126/bin/picard.jar
trimmomatic=$HOME/bin/Trimmomatic-0.32/trimmomatic-0.32.jar

# Get parameters...
MPILEUPLIST LINE
#Parameters=$(sed -n -e "$SGE_TASK_ID p" mpileupList)
sample=$( echo "$Parameters" | awk '{print $1}' )
experiment=$( echo "$Parameters" | awk '{print $2} ' )
ref=$( echo "$Parameters" | awk '{print $3}' ) ##$HOME/indexes/SL3.0/SL3.0;/sonas-hs/lippman/hpc/home/zlemmon/indexes/Pepper155_Genome.fa
snpEff_DB=$( echo "$Parameters" | awk '{print $4}' ) #### Cann_v1.55; SL3.0_ITAG3.2; SL2.50; SL2.40.21; SL3.0
proteinfile=$( echo "$Parameters" | awk '{print $5}' )
### /sonas-hs/lippman/hpc/home/zlemmon/indexes/ITAG2.4_proteins.fasta
### /sonas-hs/lippman/hpc/home/zlemmon/indexes/Pepper_1.55.gene_models.descriptions
### /sonas-hs/lippman/hpc/home/zlemmon/indexes/SL3.0/ITAG3.00/ITAG3.0_proteins.fasta
### /sonas-hs/lippman/hpc/home/zlemmon/indexes/SL3.0/ITAG3.20/ITAG3.2_proteins.fasta

# Convert sample string to array
IFS=';' read -r -a samples <<< "$sample"

####
# List of bam file locations for common parent lines etc.
####
#pimpbam=/sonas-hs/lippman/hpc/data/Spimp/SL3.0/spimp_genomic_dups.bam
#M82bam=/sonas-hs/lippman/hpc/data/Bolger2014_M82Genome_SL3.0/M82genomic_dups.bam
#utmfrenchbam=/sonas-hs/lippman/hpc/data/MicrotomGenome/SL3.0/microtom_genomic_french_dups.bam
#utmjapanesebam=/sonas-hs/lippman/hpc/data/MicrotomGenome/SL3.0/microtom_genomic_japanese_dups.bam
#s2bam=/sonas-hs/lippman/hpc/data/slike_quantification/ExD_s2_F2mapping_20170725/15-616_s2_dups.bam
#ExDbam=/sonas-hs/lippman/hpc/data/slike_quantification/ExD_s2_F2mapping_20170725/8924-2_ExD_dups.bam

# Targets if wanting to only call within gene regions.
#targetsfile=/sonas-hs/lippman/hpc/home/zlemmon/indexes/ITAG2.4_gene_regions.txt

############# Portion of script that calls SNPs among all bam files for MAPPING of interval
# Call a parallel mpileup on all samples 

# Parse samples into shorter names to see better when reporting to the log file
n=0
rm -f newsamplenames
for i in "${samples[@]}"; do
	itrunc=${i/*\//}; itrunc=${itrunc/_dups.bam/}; itrunc=${itrunc/_sorted.bam/}; sampleshort[n]=$itrunc; ((n++))
	samplestring="$samplestring, $itrunc"
	echo "$i $itrunc" >> newsamplenames
done
samplestring=${samplestring/, /}

samplenumber=$( expr ${#samples[@]} )


echo -ne "\n\n###############################################################################################\n\n"
echo "Running mpileup for $samplenumber samples: ( $samplestring )" ; date
temp_arr=($( samtools view -H "${samples[0]}" | grep "\@SQ" | sed 's/^.*SN://g' | cut -f 1 ))
common_chr_string=$(printf "%s\n" "${temp_arr[@]}" | sed -e '$!{N;s/^\(.*\).*\n\1.*$/\1\n\1/;D;}')
samtools view -H "${samples[0]}" | grep "\@SQ" | sed 's/^.*SN://g' | cut -f 1 | xargs -I LINE -n 1 -P 14 sh -c "samtools mpileup --ignore-RG -d 1000000 -t DP,AD -Q0 -Bugf "$ref" -r 'LINE' $( echo "${samples[@]}" ) | bcftools call -mvO z -o "$experiment"_'LINE'.vcf.gz"
# index all individual chromosome vcf files
echo "Indexing chromosome vcf.gz files" ; date
for file in "$experiment"_"$common_chr_string"*.vcf.gz ; do tabix -fp vcf "$file" ; done
# merge SNPs
echo "Merging SNPs called on the various chromosomes" ; date
bcftools concat -O z -o "$experiment".vcf.gz "$experiment"_"$common_chr_string"*.vcf.gz
# index final variants vcf.gz file
echo "Indexing final vcf.gz files" ; date
tabix -p vcf "$experiment".vcf.gz
echo "Done calling SNPs"
# Remove individual chromosome files.
echo -ne "Removing individual chromosome files... "
rm -f "$experiment"_"$common_chr_string"*.vcf.gz*
echo -ne 'Done!!\n\n'

# Rename samples to something more manageable than a full filename/path
echo "Renaming samples to a concise description instead of file path." ; date
bcftools reheader -s newsamplenames "$experiment".vcf.gz > "$experiment"_rehead.vcf.gz
mv -f "$experiment"_rehead.vcf.gz "$experiment".vcf.gz
bcftools index "$experiment".vcf.gz

#############
# annotate SNP calls
echo -ne "\n\n###############################################################################################\n\n"

########## Portion of script that calls effects for novel SNPs.
# Use snpEff to annotate the variants for functional effect
echo "Annotating SNPs called from experiment:$experiment." ; date
java -Xmx16g -jar $HOME/bin/snpEff/snpEff.jar ann "$snpEff_DB" -no-downstream -no-intergenic -no-upstream -stats "$experiment"_SnpEffSummary.html "$experiment".vcf.gz > "$TMPDIR"/"$experiment"_ann.vcf
bcftools view -Oz -o "$experiment"_ann.vcf.gz "$TMPDIR"/"$experiment"_ann.vcf
rm -f "$TMPDIR"/"$experiment"_ann.vcf
tabix -p vcf "$experiment"_ann.vcf.gz

# parse full file into easier to read stuff.
### input through parameter file above ###proteinfile=/sonas-hs/lippman/hpc/home/zlemmon/indexes/SL3.0/ITAG3.20/ITAG3.2_proteins.fasta
echo "Parsing "$experiment"_ann.vcf.gz into "$experiment"_ann.txt.gz" ; date
java -jar $HOME/bin/snpEff/SnpSift.jar extractFields "$experiment"_ann.vcf.gz CHROM POS REF ALT "ANN[*].GENE" "ANN[*].IMPACT" "ANN[*].EFFECT" "ANN[*].HGVS_P" | awk 'FNR==NR && />/{gsub(/>/,"",$1);gsub(/.1$/,"",$1);desc=$0;gsub(/Solyc[01][0-9]g[0-9]*.1 /,""desc);arr[$1]=desc;next}FNR==NR{FS="\t";arr[$1]=$2;next}$5 in arr{printf("%s\t%s\n",$0,arr[$5]);next}{print}' $proteinfile - | gzip > "$experiment"_ann.txt.gz

# Filter mutant variants to only include those that are homozygous variant in mutant, novel (REF in M82/pimp) and possibly on specified chromosome.
#####
#MappedChrom=SL2.50chXX # Chromosome where the putative mutation is located based on previous mapping.
# Put it in as follows...    " isRef(GEN[M82]) & ( CHROM = $MappedChrom ) & isHom(...) " ### Need to double check this later to make sure...
#echo "Filtering novel mutant SNPs for SNPs also in other haplotypes."
# Samples are: M82, pimp, e0673_mut, e0673_wt
# For e0673_mut (clb), filter so homoygous reference everything except the clb samples, which should be 1/1 and 0/1 in the mutant pools and wildtype pools, respectively.
#java -jar ~/bin/snpEff/SnpSift.jar filter " isRef(GEN[M82]) & isRef(GEN[pimp]) & ( isHom(GEN[e0673_m1_mut]) & isVariant(GEN[e0673_m1_mut]) ) & isHet(GEN[e0673_pimpF2_wt]) & isHet(GEN[e0673_m1_wt]) " "$experiment"_ann.vcf.gz | bcftools view -Oz -o e0673_mut_novelm1.vcf.gz

# Filter final VCF by a number of factors and make a file we can run SNP-index for mapping the genomic region associated with phenotype
echo -ne "\n\n###############################################################################################\n\n"
echo "Filtering for biallelic, 100bp from indel, MQ50 SNPs..." ; date
bcftools view -m2 -M2 "$experiment"_ann.vcf.gz | bcftools filter --SnpGap 100 -i ' TYPE="SNP" & MQ>=50 ' -Oz -o "$experiment"_ann_filt.vcf.gz

echo "Parsing filtered SNPs into a tab-delimited gunzip compressed text file for mapping in R." ; date
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%DP[\t%SAMPLE %GT %AD]\n' "$experiment"_ann_filt.vcf.gz | awk 'BEGIN{FS="\t"}\
        {printf("%s\t%s\t%s\t%s\t%s",$1, $2, $3, $4, $5)}
	{for(i=6;i<=NF;i++){\
		string=$i; split(string, SGD, " ");\
		sample=SGD[1] ; gt=SGD[2] ; ad=SGD[3];\
		split(ad,ad2,",");\
		printf("\t%s %s %s %s", sample, gt, ad2[1], ad2[2]);\
	}\
	printf("\n");\
}' | gzip > "$experiment"_filt.txt.gz

echo "Done with SNP calling, annotation, filtering for mapping SNPs, and filtering for novel coding effects." ; date
