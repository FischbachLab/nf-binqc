#!/usr/bin/env bash

PROJECT="${1}"
SEQKIT_DIR="${2}"
BARNNAP_DIR="${3}"
CHECKM="${4}"
GTDB="${5}"

# Get fasta stat file
for i in ${SEQKIT_DIR}/*.seqkit_stats.txt;
do
  tail -n 1 $i >> fasta_stats.txt
done

# Get rRNA counts
grep -c '>' ${BARNNAP_DIR}/*.rrna.fasta | sed -e 's/.rrna.fasta:/\t/' -e "s#${BARNNAP_DIR}/##"  >> num_rrna.txt

# Run report
binqc_report.py $PROJECT $CHECKM $GTDB fasta_stats.txt num_rrna.txt
