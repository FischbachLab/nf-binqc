# NF-BINQC

Requires that all the genomes in the `--fastas` location be in fasta format, with the same extension. By default, the pipeline assumes that they all have the extension `fna`. If this is not the case, please use the flag `--ext` to reflect the appropriate extension.

## Test

```bash
aws batch submit-job \
    --job-name nf-binqc-test \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command=FischbachLab/nf-binqc,\
"--fastas","s3://nextflow-pipelines/nf-binqc/test/data",\
"--project","00_Test",\
"--ext","fna"
```

## Actual sample command

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
