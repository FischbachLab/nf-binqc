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

bindir_ch = Channel.value("${params.fastas}")

process SEQKIT {
  tag "${id}"

  cpus 2
  memory 4.GB
  container "quay.io/biocontainers/seqkit:0.12.0--0"

  publishDir "$output_base/SeqKit/${id}"

  input:
    tuple val(id), path(assembly) from seqkit_ch

  output:
    path "${id}.seqkit_stats.txt"
    // path "*.version.txt" to versions_ch

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
//     tag "Running"

//     cpus 4
//     memory 40.GB
//     container "quay.io/biocontainers/gtdbtk:1.5.1--pyhdfd78af_0"

//     publishDir "$output_base/GTDBtk/"

//     input:
//         path(assembly_dir) from bindir_ch

//     output:
//         path "gtdbtk.${params.project}.*.summary.tsv"        
//         path "gtdbtk.${params.project}.*.classify.tree.gz"   
//         path "gtdbtk.${params.project}.*.markers_summary.tsv"
//         path "gtdbtk.${params.project}.*.msa.fasta.gz"       
//         path "gtdbtk.${params.project}.*.user_msa.fasta"     
//         path "gtdbtk.${params.project}.*.filtered.tsv"       
//         path "gtdbtk.log"                  
//         path "gtdbtk.warnings.log"         
//         path "gtdbtk.${params.project}.failed_genomes.tsv"   
//         path "*.version.txt"

//     script:
//     """
//     export GTDBTK_DATA_PATH="${params.gtdb_db}"
//     mkdir -p pplacer_tmp/${params.project}

//     gtdbtk classify_wf \\
//     --genome_dir ${assembly_dir.baseName} \\ 
//     --extension ${params.ext} \\
//     --prefix gtdb.${params.project}
//     --out_dir "\${PWD}"
//     --cpus $task.cpus \\
//     --pplacer_cpus ${params.pplacer_threads} \\
//     --scratch_dir pplacer_tmp/${params.project}

//     gzip "gtdbtk.${params.project}".*.classify.tree "gtdbtk.${params.project}".*.msa.fasta
//     gtdbtk --version | sed "s/gtdbtk: version //; s/ Copyright.*//" > gtdbtk.version.txt
//     """
// }

// TODO @sunitj: Takes and Emits a directory
process CHECKM {
    tag "Running"

    cpus 8
    memory 40.GB
    container "quay.io/biocontainers/checkm-genome:1.1.3--py_1"

    publishDir "$output_base/CheckM/"

    input:
        path(assembly_dir) from bindir_ch

    output:
        path "checkm-${params.project}.tsv"
        path "checkm-${params.project}/*"
        // path "*.version.txt" to versions_ch
// aws s3 cp ${assembly_dir}/ bins/ --exclude '*' --include '*'${params.ext}
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
      -f checkm-${params.project}.tsv \\
      ${assembly_dir.baseName} \\
      checkm-${params.project}

  ls -lhtra checkm-${params.project}/
  """
}

// Extract rRNA
process BARRNAP {
  tag "${id}"

  cpus 2
  memory 4.GB
  container "quay.io/biocontainers/barrnap:0.9--3"

  publishDir "$output_base/Barrnap/${id}"

  input:
    tuple val(id), path(assembly) from barrnap_ch

  output:
    path "${id}.rrna.gff"
    path "${id}.rrna.${params.ext}"
    // path "*.version.txt" to versions_ch

  script:
  """
    barrnap --threads $task.cpus \\
            -o ${id}.rrna.${params.ext} \\
            $assembly > ${id}.rrna.gff
    barrnap --version > barnap.version.txt
  """
}


// Collect versions; concatenate files; remove duplicates; save