#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// If the user uses the --help flag, print the help text below
params.help = false

// Function which prints help message text
def helpMessage() {
  log.info"""
    QC for MAGs

    Required Arguments:
      --fastas   path       Location where 1 or more fasta files are stored.
      --project  value      Folder to place analysis outputs (default: )
      --seedfile path       A csv file with headers in the format "Name,Fasta",

    Options
      --ext       value     Extension of the fasta files in the fastas (default:)
      --outdir    path      Base directory for output files (default: )
      --tree  <true|false>  Optional whether a tree is built or not (default:false)
    """.stripIndent()
}

log.info"""Starting""".stripIndent()

// Show help message if the user specifies the --help flag at runtime
if (params.help) {
  // Invoke the function above which prints the help message
  helpMessage()
  // Exit out and do not run anything else
  exit 0
}

// Check input options
if (params.fastas && params.seedfile){
   exit 1, "Input fasta files must be defined using either '--fastas' or '--seedfile' parameter. Please choose one way"
}

def outputBase = "${params.outdir}/${params.project}"

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
        publishDir "$outputBase"

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


// SEQKIT
process SEQKIT {
  tag "${id}"

  container params.docker_container_seqkit

  publishDir "$outputBase/01_SeqKit/${id}"

  input:
    tuple val(id), path(assembly)

  output:
    path "${id}.seqkit_stats.txt", emit: seqkit_out_ch
    path "*.version.txt"

  script:
  """
  seqkit \\
    stats \\
    -T \\
    -j $task.cpus \\
    $assembly 1> ${id}.seqkit_stats.txt
  seqkit version > seqkit.version.txt
  """
}

// // Extract rRNA
process BARRNAP {
  tag "${id}"

  container params.docker_container_barrnap

  publishDir "$outputBase/02_Barrnap/${id}"

  input:
    tuple val(id), path(assembly)

  output:
    path "${id}.rrna.gff"
    path "${id}.rrna.fasta", emit: barrnap_out_ch
    path "*.version.txt"

  script:
  """
    barrnap --threads $task.cpus \\
            -o ${id}.rrna.fasta \\
            $assembly > ${id}.rrna.gff
    barrnap -v > barnap.version.txt
  """
}

// CheckM
process CHECKM {
  tag "${params.project}"

  container params.docker_container_checkm

  publishDir "$outputBase/03_CheckM"

  input:
    path 'assembly_dir'

  output:
    path "checkm-lineage.tsv"
    path "checkm-qa.tsv", emit: checkm_out_ch
    path "*.version.txt"


  script:
  """
  checkm data setRoot ${params.checkm_db}
  checkm \\
    lineage_wf \\
      --nt \\
      --ali \\
      --tab_table \\
      -t $task.cpus \\
      -x ${params.ext} \\
      --pplacer_threads ${params.pplacer_threads} \\
      -f checkm-lineage.tsv \\
      assembly_dir \\
      checkm

  checkm \\
    qa \\
      -o 2 \\
      --tab_table checkm/lineage.ms \\
      checkm 1> checkm-qa.tsv 2> checkm-qa.log

  echo "CheckM" > checkm.version.txt
  checkm &> checkm.version.txt
  """
// tar cvzf checkm-${params.project}.tar.gz checkm-${params.project}
}


// GTDBtk
process GTDBTK {
    tag "${params.project}"

    container params.docker_container_gtdbtk

    publishDir "$outputBase/04_GTDBtk"

    input:
      path 'assembly_dir'

    output:
      // path "*.summary.tsv"
      // path "*.classify.tree"
      // path "*.markers_summary.tsv"
      // path "*.msa.fasta"
      // path "*.user_msa.fasta"
      // path "*.filtered.tsv"
      // path 'gtdbtk.log'
      // path 'gtdbtk.warnings.log'
      // path "*.failed_genomes.tsv"
      path "*.version.txt"
      path "gtdbtk-results/*"
      path  "gtdbtk-results/gtdb.${params.project}.bac120.summary.tsv", emit: gtdb_out_ch

    script:
    // def pplacer_threads = half(${task.cpus})
    """
    export GTDBTK_DATA_PATH="${params.gtdb_db}"
    mkdir -p pplacer_tmp/${params.project}

    gtdbtk classify_wf \\
      --genome_dir assembly_dir \\
      --extension ${params.ext} \\
      --prefix gtdb.${params.project} \\
      --out_dir gtdbtk-results \\
      --cpus $task.cpus \\
      --pplacer_cpus ${params.pplacer_threads} \\
      --skip_ani_screen

    echo "GTDBtk" > gtdbtk.version.txt
    gtdbtk --version | sed "s/gtdbtk: version //; s/ Copyright.*//" >> gtdbtk.version.txt
    echo "GTDBTK_DATA_PATH = \${GTDBTK_DATA_PATH}" >> gtdbtk.version.txt
    """
    // --scratch_dir pplacer_tmp/${params.project}
    // gzip "gtdbtk.${params.project}".*.classify.tree
    // gzip "gtdbtk.${params.project}".*.msa.fasta
}

// // Collect versions; concatenate files; remove duplicates; save
// versions_ch
//     .collectFile(name: out)
// Report
process REPORT {
    tag "${params.project}"

    container params.docker_container_report

    publishDir "$outputBase/05_REPORT"
    publishDir "s3://genomics-workflow-core/aws-miti-straindb-us-west-2/aws_glue/assembly_qc/"

input:
      path "seqkit_dir/*"
      path "barrnap_rrna_dir/*"
      path checkm_qa
      path summary
output:
      path "${params.project}.report.csv"
script:
"""
     bash qc_report_wrapper.sh "${params.project}" seqkit_dir barrnap_rrna_dir $checkm_qa $summary
"""

}


process kSNP4_tree {
    tag "${params.project}"

    errorStrategy 'ignore'
    container params.docker_container_tree

    publishDir "$outputBase/06_TREE"

input:
      path "genomes"

output:
      file "kSNP4_output/tree.SNPs_all.ML.tre"
script:

"""
  MakeKSNP4infile -indir genomes -outfile genomes.list
  kSNP4 -in genomes.list -k 21 -ML -outdir kSNP4_output
"""

}

/*
 Write the tree in pdf using figtree
*/
 process write_tree_in_pdf {
     tag "${params.project}"

     container params.docker_container_figtree

     publishDir "$outputBase/06_TREE", mode:'copy', pattern: '*.pdf'

     input:
     file tree

     output:
     file "*.pdf"

     script:
     """
     figtree -graphic PDF "$tree" "${params.project}".bestTree.pdf
     """
 }

workflow {

  if(params.seedfile){

   Channel.fromPath(params.seedfile) | save_seedfile

   seqkit_barrnap_ch = Channel
            .fromPath(params.seedfile)
            .ifEmpty { exit 1, "Cannot find any seed file matching: ${params.seedfile}." }
            .splitCsv(header: ['name', 'file'], sep: ',')
            .map{ row -> tuple(row.name, file(row.file))}

    // Save a copy of input fasta files
    seqkit_barrnap_ch | copy_fastas

    copy_fastas.out.copy_fastas_ch | dummy
    gtdb_checkm_bindir_ch = copy_fastas.out.copy_fastas_ch
    // tmp fix
    gtdb_checkm_bindir_ch = Channel.fromPath("$outputBase/00_Fasta")
    //gtdb_checkm_bindir_ch = Channel
    //       .fromPath(params.seedfile)
    //       .splitCsv(header: ['name', 'file'], sep: ',')
    //       .map{ row -> file(row.file)}

  } else {

     // Show help message if the user specifies a fasta file but not makedb or db
     if (params.fastas == "") {
      // Invoke the function above which prints the help message
      helpMessage()
      // Exit out and do not run anything else
      exit 0
    }
      /*
       * Defines the pipeline inputs parameters (giving a default value for each for them)
       * Each of the following parameters can be specified as command line options
       */

      def fnaGlob = "${params.fastas}/*.${params.ext}"
      log.info"""Searching for file at this location: $fnaGlob""".stripIndent()

      // Exit if the file is empty
      Channel
        .fromPath(fnaGlob)
        .ifEmpty { exit 1, 'Cannot find matching fasta file' }

      /*
        Since the processes are not dependent upon each other, duplicate the channel
        such that there is one channel for each process
       */
      seqkit_barrnap_ch = Channel
        .fromPath(fnaGlob)
        .map { it -> tuple(it.baseName, file(it)) }

      // Pass the entire directory
      gtdb_checkm_bindir_ch = Channel
        .fromPath(params.fastas)

  }


        seqkit_barrnap_ch |  SEQKIT
        seqkit_barrnap_ch | BARRNAP
        gtdb_checkm_bindir_ch | CHECKM
        gtdb_checkm_bindir_ch | GTDBTK
        REPORT( SEQKIT.out.seqkit_out_ch.toSortedList(), BARRNAP.out.barrnap_out_ch.toSortedList(), CHECKM.out.checkm_out_ch, GTDBTK.out.gtdb_out_ch)


    if (params.tree){
        gtdb_checkm_bindir_ch | kSNP4_tree
        kSNP4_tree.out | write_tree_in_pdf

    }
}
