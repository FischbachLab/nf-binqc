#!/usr/bin/env nextflow
nextflow.enable.dsl=1
// If the user uses the --help flag, print the help text below
params.help = false

// Function which prints help message text
def helpMessage() {
  log.info"""
    QC for MAGs

    Required Arguments:
      --fastas              Location where 1 or more fasta files are stored.
      --project             Folder to place analysis outputs (default: )

    Options
      --ext                 Extension of the fasta files in the fastas (default:)
      --outdir              Base directory for output files (default: )
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

// // Show help message if the user specifies a fasta file but not makedb or db
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

def outputBase = "${params.outdir}/${params.project}"

/*
  Since the processes are not dependent upon each other, duplicate the channel
  such that there is one channel for each process
 */
Channel
  .fromPath(fnaGlob)
  .map { it -> tuple(it.baseName, file(it)) }
  .into { seqkit_ch; barrnap_ch }

// Pass the entire directory
Channel
  .fromPath( "${params.fastas}", type: 'dir')
  // .set{checkm_bindir_ch}
  // .set{gtdb_bindir_ch}
  .into { gtdb_bindir_ch; checkm_bindir_ch}

// checkm_bindir_ch.view()
// gtdb_bindir_ch.view()

// SEQKIT
process SEQKIT {
  tag "${id}"

  container params.docker_container_seqkit

  publishDir "$outputBase/01_Stats/${id}"

  input:
    tuple val(id), path(assembly) from seqkit_ch

  output:
    path "${id}.seqkit_stats.txt"
    path "*.version.txt"
    path "${id}.sha256"

  script:
  """
  seqkit \\
    stats \\
    -T \\
    -j $task.cpus \\
    $assembly 1> ${id}.seqkit_stats.txt
  seqkit version > seqkit.version.txt

  sha256sum $assembly &> ${id}.sha256
  """
}

// // Extract rRNA
process BARRNAP {
  tag "${id}"

  container params.docker_container_barrnap

  publishDir "$outputBase/02_RRNA/${id}"

  input:
    tuple val(id), path(assembly) from barrnap_ch

  output:
    path "${id}.rrna.gff"
    path "${id}.rrna.${params.ext}"
    path "*.version.txt"

  script:
  """
    barrnap --threads $task.cpus \\
            -o ${id}.rrna.${params.ext} \\
            $assembly > ${id}.rrna.gff
    barrnap -v > barnap.version.txt
  """
}

// CheckM
process CHECKM {
  tag "${params.project}"

  container params.docker_container_checkm

  publishDir "$outputBase/03_Contamination"

  input:
    path(assembly_dir) from checkm_bindir_ch

  output:
    path "checkm-lineage.tsv"
    path "checkm-qa.tsv"
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
      ${assembly_dir} \\
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

    publishDir "$outputBase/04_Classification"

    input:
      path(assembly_dir) from gtdb_bindir_ch

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


    script:
    // def pplacer_threads = half(${task.cpus})
    """
    export GTDBTK_DATA_PATH="${params.gtdb_db}"
    mkdir -p pplacer_tmp/${params.project}
  
    gtdbtk classify_wf \\
      --genome_dir ${assembly_dir} \\
      --extension ${params.ext} \\
      --prefix gtdb.${params.project} \\
      --out_dir gtdbtk-results \\
      --cpus $task.cpus \\
      --pplacer_cpus ${params.pplacer_threads}

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

// TODO: @sunitj
// replicate bin/qc_wrapper.sh as a Nextflow pipeline 