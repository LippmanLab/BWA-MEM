#!/bin/bash

if [ $# -ne 3 ] ; then echo -ne "Usage: ./prepScripts.sh ProjectName ParameterFile mpileupList\n\tProjectName: unique string identifying project\n\tParameterFile: tab-delimited file carrying parameters for bwa mem alignment\n\tmpileupList: single line tab-delim file with bam alignments and experiment string\n\n"; exit 0; fi

# Make input easier to deal with...
Proj=$1
PF=$2
mL=$3

# Make a log directory path and store to pass to next files...
echo "Making log dir..."
mkdir -p log
logdir=$PWD/log

# Count SampleNumber
echo "Counting samples for bwa mem..."
SampleNumber=( `wc -l $PF | cut -f1` )

# Adjust log line and parallel line in the two scripts to have the appropriate number of parallel processes etc...
echo "Making align_$Proj.sh and mpileup_$Proj.sh"
awk -v SampleNumber="$SampleNumber" -v logdir="$logdir" -v ParameterFile="$PF" '{if($0 ~ /^SAMPLENUMBER LINE$/) {print "#$ -t 1-" SampleNumber} else if ($0 ~ /^LOGDIR LINE$/) {print "#$ -o " logdir} else if ($0 ~ /^PARAMETERFILE LINE$/) {print "Parameters=$(sed -n -e \"$SGE_TASK_ID p\" " ParameterFile ")"} else {print}}' align.sh > align_$Proj.sh
awk -v logdir="$logdir" -v mpileupList="$mL" '{if($0 ~ /^LOGDIR LINE$/) {print "#$ -o " logdir} else if ($0 ~ /^MPILEUPLIST LINE$/) {print "Parameters=$(sed -n -e \"$SGE_TASK_ID p\" " mpileupList ")"} else {print}}' mpileup.sh > mpileup_$Proj.sh
