#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include {save_seedfile; dummy; copy_fastas } from './modules/house_keeping'
include { SEQKIT; BARRNAP; CHECKM; GTDBTK; REPORT; GUNC; kSNP4_tree; write_tree_in_pdf} from './modules/bin_qc'

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



workflow {

  if(params.seedfile){

   Channel.fromPath(params.seedfile) | save_seedfile

   seqkit_barrnap_ch = Channel
            .fromPath(params.seedfile)
            .ifEmpty { exit 1, "Cannot find any seed file matching: ${params.seedfile}." }
            .splitCsv(header: ['name', 'file'], sep: ',')
            .map{ row -> tuple(row.name, row.file)}

    // Save a copy of input fasta files
    seqkit_barrnap_ch | copy_fastas

    copy_fastas.out.copy_fastas_ch.toSortedList() | dummy

    seqkit_barrnap_ch = copy_fastas.out.copy_fastas_ch
        .map { it -> tuple(it.baseName, file(it)) }

    gtdb_checkm_bindir_ch = copy_fastas.out.copy_fastas_ch.toSortedList()

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
        gtdb_checkm_bindir_ch | GUNC
        REPORT(SEQKIT.out.seqkit_out_ch.toSortedList(), BARRNAP.out.barrnap_out_ch.toSortedList(), CHECKM.out.checkm_out_ch, GTDBTK.out.gtdb_out_ch, GUNC.out.gunc_out_ch)
        
    if (params.tree){
        gtdb_checkm_bindir_ch | kSNP4_tree
        kSNP4_tree.out | write_tree_in_pdf

    }
}
