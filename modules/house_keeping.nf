// save a copy of the seedfile
    process save_seedfile {

        container params.docker_container_report
        publishDir "s3://genomics-workflow-core/aws-miti-straindb-us-west-2/aws_glue/qc_seedfiles/"

        input:
        path seed
        output:
        path "${seed}"

       script:
       """
        ls $seed
       """
    }

    process dummy {

        container params.docker_container_seqkit

        input:
        path ('in_dir')

        output:

        script:
        """
        ls -l in_dir
        """
    }

    process copy_fastas {
        tag "$id"

        //storeDir './files'
        //stageOutMode 'copy'
        publishDir "${params.outdir}/${params.project}"

        container params.docker_container_seqkit

        input:
        tuple val(id), path(assembly)

        output:
        //path "00_Fasta/${id}.${params.ext}"
        path "00_Fasta", emit: copy_fastas_ch
        //tuple val(id), path("${id}.${params.ext}")

        script:
        """
        mkdir 00_Fasta
        cp $assembly 00_Fasta/${id}.${params.ext}
        """
    }