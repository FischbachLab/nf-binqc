NF-BIN-QC
====================

```bash
aws batch submit-job \
    --profile maf \
    --job-name nf-binqc \
    --job-queue default-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command=s3://nextflow-pipelines/nf-annotate,\
"--fastas","s3://dev-scratch/ReferenceDBs/NinjaMap/Index/12Com_20210528/fasta",\
"--project","12Com_20210528"
```
