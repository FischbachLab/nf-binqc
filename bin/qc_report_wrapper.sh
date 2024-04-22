#!/usr/bin/env bash

PROJECT="${1}"
SEQKIT_DIR="${2}"
BARNNAP_DIR="${3}"
CHECKM="${4}"
GTDB="${5}"
GUNC="${6}"

# Get fasta stat file
for i in ${SEQKIT_DIR}/*.seqkit_stats.txt;
do
  tail -n 1 $i >> fasta_stats.txt
done

# Get rRNA counts
line=`ls  ${BARNNAP_DIR}/*.rrna.fasta | wc -l`
if [ $line -gt 1 ]
then
     grep -c '>' ${BARNNAP_DIR}/*.rrna.fasta | sed -e 's/.rrna.fasta:/\t/' -e "s#${BARNNAP_DIR}/##"  >> num_rrna.txt
# handle a single input file
else
     file=`ls ${BARNNAP_DIR}/*.rrna.fasta`
     trimmed_file1=${file##*/}
     trimmed_file2=${trimmed_file1%%.*}
     num=`grep -c '>' ${BARNNAP_DIR}/*.rrna.fasta`
     echo -e $trimmed_file2"\t"$num > num_rrna.txt
fi


# Run report
binqc_report.py $PROJECT $CHECKM $GTDB fasta_stats.txt num_rrna.txt $GUNC
