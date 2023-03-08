NF-BINQC
====================

```bash
aws batch submit-job \
    --job-name nf-binqc-1019-1 \
    --job-queue default-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command=s3://nextflow-pipelines/nf-binqc,\
"--fastas","s3://nextflow-pipelines/nf-binqc/test/data"
```

```bash
aws batch submit-job \
    --job-name nf-binqc-SCv2_4_20210212 \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command=s3://nextflow-pipelines/nf-binqc,\
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
## Note that all assembled genome files are saved in one location and supplied to the --fastas option
### One step submission to generate report. No local post-processing needed any longer.

```bash
aws batch submit-job \
    --job-name nf-binqc-MCB \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command="FischbachLab/nf-binqc, \
"--ext", "fasta", \
"--fastas","s3://genomics-workflow-core/Results/BinQC/MITI-MCB/20221018/fasta", \
"--project","20221018_207_v2", \
"--outdir","s3://genomics-workflow-core/Results/BinQC/MITI-MCB" "
```


### The QC report is in the 05_REPORT directory.
```bash
s3://genomics-workflow-core/Results/BinQC/MITI-MCB/20221018_207_v2/05_REPORT/20221018_207_v2.report.csv
```
