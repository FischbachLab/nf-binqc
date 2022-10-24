#!/usr/bin/env bash

# generate QC report after running nf-binqc 

DIR=${1:?"Specify a dir path as argv[1]"}
DB=${2:?"Specify a DB name as argv[2]"}

mkdir -p $DIR

echo $DIR

cd $DIR

#Download qc results
aws s3 sync s3://genomics-workflow-core/Results/BinQC/MITI-MCB/$DB/ qc/

if [ -d qc/stats ]; then
   rm -r qc/stats/
fi

mkdir -p qc/stats/

# get fasta stat file
ls qc/01_SeqKit/*/*.seqkit_stats.txt | parallel "tail -n 1 {} >> qc/stats/fasta_stats.txt"

# Keep Genes and Proteins
#mkdir -p qc/orfs/{genes,proteins}; \
#find qc/04_GTDBtk/ -name '*_protein.fna' | parallel "cp {} qc/orfs/genes/"; \
#find qc/04_GTDBtk/ -name '*_protein.faa' | parallel "cp {} qc/orfs/proteins/"; \
#grep -c '>' qc/orfs/genes/*_protein.fna | sed -e 's/_protein.fna:/\t/' -e "s#qc/orfs/genes/##" > qc/orfs/num_genes.txt

#Get rRNA counts
mkdir -p qc/stats/rrna; \
find qc/02_Barrnap/ -name '*.rrna.fasta' | parallel "cp {} qc/stats/rrna/"; \
grep -c '>' qc/stats/rrna/*.rrna.fasta | sed -e 's/.rrna.fasta:/\t/' -e "s#qc/stats/rrna/##"  > qc/stats/num_rrna.txt 


# Run report
#python generate_ninjamap_db_report.py input_dir_binqc db_name
python ~/scripts/generate_binqc_report.py $DIR $DB
