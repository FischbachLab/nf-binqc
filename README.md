NF-BINQC
====================

```bash
aws batch submit-job \
    --job-name nf-binqc-1019-1 \
    --job-queue default-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command=FischbachLab/nf-binqc,\
"--fastas","s3://nextflow-pipelines/nf-binqc/test/data"
```

```bash
aws batch submit-job \
    --job-name nf-binqc-SCv2_4_20210212 \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command=FischbachLab/nf-binqc,\
"--fastas","s3://maf-versioned/ninjamap/Index/SCv2_4_20210212/fasta",\
"--project","SCv2_4_20210212"
```

or locally

```bash
nextflow run . --fastas test/data --project 00_TEST --outdir test/result/
```

## Update

Note that updating the code here will *not* update the pipeline automatically.

```{bash}
cd nf-binqc
aws s3 sync . s3://nextflow-pipelines/nf-binqc --exclude ".git/*" --exclude "test/result/*" --delete --profile maf
```

# MITI MCB strain QC example batch job
### One step submission to generate report. No local post-processing needed any longer.
### Input fasta files must be defined using either '--fastas' or '--seedfile' parameter.
### Note that all assembled genome files are saved in one location and supplied to the --fastas option or supply a seedfile defining the genome name and s3 path in a csv file without headers. Below is an example seedfile.  

```
SH0001342-00095,s3://genomics-workflow-core/Results/HybridAssembly/MITI-MCB/SH0001342-00095/UNICYCLER/assembly.fasta
SH0001372-00039,s3://genomics-workflow-core/Results/HybridAssembly/MITI-MCB/SH0001372-00039/UNICYCLER/assembly.fasta
```

#### The batch submission example using --fastas option

```{bash}
aws batch submit-job \
    --job-name nf-binqc-MCB \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command="FischbachLab/nf-binqc, \
"--ext", "fasta", \
"--fastas","s3://genomics-workflow-core/Results/BinQC/MITI-MCB/20221018/fasta", \
"--project","20221018_207_v2", \
"--tree","tree", \
"--outdir","s3://genomics-workflow-core/Results/BinQC/MITI-MCB" "
```

#### The batch submission example using --seedfile option

```{bash}
aws batch submit-job \
    --job-name nf-binqc-MCB \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command="FischbachLab/nf-binqc, \
"--ext", "fasta", \
"--seedfile", "s3://genomics-workflow-core/Results/BinQC/TEST/20221018_207_v2.seedfile.csv"
"--project","20221018_207_v2", \
"--tree","true", \
"--outdir","s3://genomics-workflow-core/Results/BinQC/MITI-MCB" "
```

#### The structure of the output directory
```
00_Fasta/
01_SeqKit/
02_Barrnap/
03_CheckM/
04_GTDBtk/
05_REPORT/
06_TREE/
```
### Example: The QC report and pdf tree files are saved in the 05_REPORT and 06_TREE directories, respectively.
```
s3://genomics-workflow-core/Results/BinQC/MITI-MCB/20221018_207_v2/05_REPORT/20221018_207_v2.report.csv
s3://genomics-workflow-core/Results/BinQC/MITI-MCB/20221018_207_v2/06_TREE/20221018_207_v2.bestTree.pdf
```

### Reference: Generating the Whole Genome Trees
```
[kSNP4.1](https://pubmed.ncbi.nlm.nih.gov/37948764/) is a k-mer based SNP calling tool that can be used to generate whole genome trees.
```
