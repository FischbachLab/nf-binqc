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
  .fromPath("${params.fastas}/*.${params.ext}", type : 'file', checkIfExists: true)
  .ifEmpty { exit 1, "Cannot find matching fasta file" }

def output_base = "${params.outdir}/${params.project}"
/* 
  Since the processes are not dependent upon each other, duplicate the channel 
  such that there is one channel for each process
 */
Channel
    .fromPath(fna_glob, type : 'file')
    .map { file -> [file.baseName, file] }
    .into {seqkit_ch; barrnap_ch}

// Channel
//     .fromPath("${params.fastas}", type:'dir')
//     into{checkm_ch}

checkm_ch = Channel.value("${params.fastas}")


// Channel
//     .fromPath(params.fastas, type : 'dir')
//     // .into {seqkit_ch; gtdbtk_ch; checkm_ch; barrnap_ch; pgap_ch}
//     .into {gtdbtk_ch; checkm_ch}

process SEQKIT {
  tag "Genome Stats: ${id}"

  cpus 2
  memory 4.GB
  container "quay.io/biocontainers/seqkit:0.12.0--0"

  publishDir "$output_base/SeqKit/${id}"

  input:
    tuple val(id), path(assembly) from seqkit_ch

  output:
    path "${id}.seqkit_stats.txt"
    // path '*.version.txt' to versions_ch

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

// TODO @sunitj: Takes and Emits a directory
// process GTDBTK {
//   tag "Genome Stats: ${assembly.baseName}"

//   cpus 2
//   memory 8.GB
//   // TODO @sunitj: Will need to create a new docker with a different entrypoint
//   container "quay.io/biocontainers/seqkit:0.12.0--0"

//   publishDir "$output_base/GTDBtk/"

// input:
//     tuple val(meta), path("bins/*")
//     tuple val(db_name), path("database/*")

// output:
// path "gtdbtk.${meta.assembler}-${meta.id}.*.summary.tsv"        , emit: summary
// path "gtdbtk.${meta.assembler}-${meta.id}.*.classify.tree.gz"   , emit: tree
// path "gtdbtk.${meta.assembler}-${meta.id}.*.markers_summary.tsv", emit: markers
// path "gtdbtk.${meta.assembler}-${meta.id}.*.msa.fasta.gz"       , emit: msa
// path "gtdbtk.${meta.assembler}-${meta.id}.*.user_msa.fasta"     , emit: user_msa
// path "gtdbtk.${meta.assembler}-${meta.id}.*.filtered.tsv"       , emit: filtered
// path "gtdbtk.${meta.assembler}-${meta.id}.log"                  , emit: log
// path "gtdbtk.${meta.assembler}-${meta.id}.warnings.log"         , emit: warnings
// path "gtdbtk.${meta.assembler}-${meta.id}.failed_genomes.tsv"   , emit: failed
// path '*.version.txt'                                            , emit: version

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
process CHECKM {
    tag "Running CheckM"

    cpus 4
    memory 40.GB
    container "quay.io/biocontainers/checkm-genome:1.1.3--py_1"

    publishDir "$output_base/CheckM/"

    input:
        path(assembly_dir) from checkm_ch

    output:
        path "checkm.tsv"
        path "checkm_output.tar.gz"
        // path '*.version.txt' to versions_ch
// aws s3 cp ${assembly_dir}/ bins/ --exclude '*' --include '*'${params.ext}
  script:
  """
  ls -lhtra
  
  checkm data setRoot ${params.checkm_db}
  checkm \\
    lineage_wf \\
      --nt \\
      --ali \\
      --tab_table \\
      -t $task.cpus \\
      -x ${params.ext} \\
      --pplacer_threads ${params.pplacer_threads} \\
      -f checkm.tsv \\
      . \\
      checkm_output

  tar czf checkm_output.tar.gz checkm_output
  """
//   checkm | head -n 2 | sed "/^$/d" | cut -d : -f 4 | sed "s/^ //" > checkm.version.txt
}

// Extract rRNA
process BARRNAP {
  tag "Extracting 16S for: ${id}"

  cpus 2
  memory 4.GB
  container "quay.io/biocontainers/barrnap:0.9--3"

  publishDir "$output_base/Barrnap/${id}"

  input:
    tuple val(id), path(assembly) from barrnap_ch

  output:
    path "${id}.rrna.gff"
    path "${id}.rrna.${params.ext}"
    // path '*.version.txt' to versions_ch

  script:
  """
    barrnap --threads $task.cpus \\
            -o ${id}.rrna.${params.ext} \\
            $assembly > ${id}.rrna.gff
    barrnap --version > barnap.version.txt
  """
}


// Collect versions; concatenate files; remove duplicates; save