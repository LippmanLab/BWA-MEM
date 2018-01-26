# BWA-MEM
Pipeline implemented by Zachary H Lemmon <zlemmon@cshl.edu> Oct 30, 2017. This set of scripts can be used to set up an SGE qsub for parallel processing of Illumina WGS data to a given reference.

Two parameter files needed.
1) ParameterFile = multi-line file with semi-colon separated list of "Read1a;Read1b	Read2a;Read2b	Sample_a;Sample_b	Library	BWAIndex_Path"
2) mpileupList = semi-colon separated single line file listing all bam alignment files wanting to do the mpileup with and a library/experiment name "bam1;bam2;bam3;bam4	experiment"
