NF-BINQC
====================

```bash
aws batch submit-job \
    --profile maf \
    --job-name nf-binqc-0827-1 \
    --job-queue default-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command=s3://nextflow-pipelines/nf-binqc,\
"--fastas","s3://dev-scratch/fasta/genomes"
```

or locally

```bash
nextflow run . --fastas s3://dev-scratch/fasta/genomes --project 12Com_20210528
```
