#!/usr/bin/env nextflow

// If the user uses the --help flag, print the help text below
params.help = false

// Function which prints help message text
def helpMessage() {
    log.info"""
    QC for MAGs
    
    Required Arguments:
      --fastas              Location where 1 or more fasta files are stored.
      --project             Folder to place analysis outputs (default: ${params.project})
    
    Options
      --ext                 Extension of the fasta files in the fastas (default: ${params.ext})
      --outdir              Base directory for output files (default: ${params.outdir})
    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

// // Show help message if the user specifies a fasta file but not makedb or db
if ((params.fastas  == null)){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

/*
 * Defines the pipeline inputs parameters (giving a default value for each for them) 
 * Each of the following parameters can be specified as command line options
 */

def fna_glob = "${params.fastas}" + "/*." + "${params.ext}"
log.info"""Searching for file at this location: $fna_glob""".stripIndent()

// Exit if the file is empty
Channel
  .fromPath(fna_glob, type : 'file', checkIfExists: true)
  .ifEmpty { exit 1, "Cannot find matching fasta file" }

def output_base = "${params.outdir}/${params.project}"
/* 
  Since the processes are not dependent upon each other, duplicate the channel 
  such that there is one channel for each process
 */
Channel
    .fromPath(fna_glob, type : 'file')
    // .into {seqkit_ch; gtdbtk_ch; checkm_ch; barrnap_ch; pgap_ch}
    .into {seqkit_ch; barrnap_ch}

// Channel
//     .fromPath(params.fastas, type : 'dir')
//     // .into {seqkit_ch; gtdbtk_ch; checkm_ch; barrnap_ch; pgap_ch}
//     .into {gtdbtk_ch; checkm_ch}

process SEQKIT {
  tag "Genome Stats: ${assembly.baseName}"

  cpus 2
  memory 8.GB
  container "quay.io/biocontainers/seqkit:0.12.0--0"

  publishDir "$output_base/SeqKit/${assembly.baseName}"
  

  input:
    path assembly from seqkit_ch

  output:
    path "${assembly.baseName}.seqkit_stats.txt"

  script:
  """
  seqkit \\
    stats \\
    -T \\
    -j $task.cpus \\
    $assembly 1> ${assembly.baseName}.seqkit_stats.txt
  """
}

// TODO @sunitj: Takes and Emits a directory
// process GTDBTK {
//   tag "Genome Stats: ${assembly.baseName}"

//   cpus 2
//   memory 8.GB
//   // TODO @sunitj: Will need to create a new docker with a different entrypoint
//   container "quay.io/biocontainers/seqkit:0.12.0--0"

//   publishDir "$output_base/GTDBtk/"

//   input:
//     path bin_dir from gtdbtk_ch

//   output:
//     path "${sample_name}.seqkit_stats.txt"

//   script:
//   """
//   classify_wf \\
//     --cpus $task.cpus \\
//     --extension ${params.ext} \\
//     --genome_dir $bin_dir \\ ## TODO: This is run on an entire directory!!!
//     --out_dir  ## TODO: This outputs an entire directory!!!
//   """
// }

// TODO @sunitj: Takes and Emits a directory
// process CHECKM {
//   publishDir "$output_base/${assembly.baseName}/CheckM", mode: 'move', overwrite: false

//   input:
//     path assembly from checkm_ch

//   output:
//     path "${sample_name}.seqkit_stats.txt"

//   script:
//   """
//   checkm \\
//     lineage_wf \\
//       --nt \\
//       --ali \\
//       --tab_table \\
//       -t ${THREADS} \\
//       -x ${BIN_FASTA_EXT} \\
//       --pplacer_threads ${PP_THREADS} \\
//       -f ${CHECKM_OUTPUT_FILE} \\
//       ${BIN_DIR} \\ ## TODO: This is run on an entire directory!!! 
//       ${CHECKM_OUTPUT_DIR} ## TODO: This outputs an entire directory!!!
//   """
// }

// Extract rRNA
process BARRNAP {
  tag "Extracting 16S for: ${assembly.baseName}"

  cpus 2
  memory 8.GB
  container "quay.io/biocontainers/barrnap:0.9--3"

  publishDir "$output_base/Barrnap/${assembly.baseName}"

  input:
    path assembly from barrnap_ch

  output:
    path "${assembly.baseName}.rrna.gff"
    path "${assembly.baseName}.rrna.${params.ext}"


  script:
  """
  barrnap --threads $task.cpus \\
            -o ${assembly.baseName}.rrna.${params.ext} \\
            $assembly > ${assembly.baseName}.rrna.gff
  """
}
