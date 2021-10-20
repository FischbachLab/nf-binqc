#!/bin/bash

#sample command
#/home/ec2-user/efs/docker/Xmeng/qc_wrapper.sh /Path/to/genomes db_name

Genomes=${1:?"Specify a dir path as argv[1]"}
DB=${2:?"Specify a DB name as argv[2]"}
echo "Fasta path: $Genomes"
echo "DB name: $DB"

# build a seedfile
cd $Genomes
echo "Strain_Name,FTP_link" > ../$DB.seedfile.txt && \
for i in *; do echo "${i%.fna},$Genomes/$i" >> ../$DB.seedfile.txt; done

cd -

# Run generate DB
python /home/ec2-user/efs/docker/Xmeng/create_ninjamap_db.py -list $DB.seedfile.txt -db $DB

# Run Genome QC
# bash -x ../../run_binqc.sh OUTPUT_DIR_NAME  INPUT_FASTA_DIR FASTA_EXTENSION &> run_binqc.log

bash -x /home/ec2-user/efs/docker/Xmeng/run_binqc.sh qc $DB/fasta fna &> run_binqc.log 

# Keep Genes and Proteins
mkdir -p qc/orfs/{genes,proteins}; \
find qc/GTDBtk/ -name '*_protein.fna' | parallel "cp {} qc/orfs/genes/"; \
find qc/GTDBtk/ -name '*_protein.faa' | parallel "cp {} qc/orfs/proteins/"; \
grep -c '>' qc/orfs/genes/*_protein.fna | sed -e 's/_protein.fna:/\t/' -e "s#qc/orfs/genes/##" > qc/orfs/num_genes.txt

# Run report
#python generate_ninjamap_db_report.py input_dir_binqc db_name
python /home/ec2-user/efs/docker/Xmeng/generate_report.py . $DB
