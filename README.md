NF-BINQC
====================

## Test 

```bash
aws batch submit-job \
    --profile maf \
    --job-name nf-binqc-1019-1 \
    --job-queue default-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command=s3://nextflow-pipelines/nf-binqc,\
"--fastas","s3://nextflow-pipelines/nf-binqc/test/data"
```

## Actual sample command

```bash
aws batch submit-job \
    --profile maf \
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
