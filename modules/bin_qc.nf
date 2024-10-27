// SEQKIT
process SEQKIT {
  tag "${id}"

  container params.docker_container_seqkit

  publishDir "${params.outdir}/${params.project}/01_SeqKit/${id}"

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
  sed -i "s/^${id}.${params.ext}/${id}/" ${id}.seqkit_stats.txt
  """
}

// Extract rRNA
process BARRNAP {
  tag "${id}"

  container params.docker_container_barrnap

  publishDir "${params.outdir}/${params.project}/02_Barrnap/${id}"

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

  publishDir "${params.outdir}/${params.project}/03_CheckM"

  input:
    path 'assembly_dir/*'

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

    publishDir "${params.outdir}/${params.project}/04_GTDBtk"

    input:
      path 'assembly_dir/*'

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

    publishDir "${params.outdir}/${params.project}/05_REPORT"
    publishDir "s3://genomics-workflow-core/aws-miti-straindb-us-west-2/aws_glue/assembly_qc/"

input:
      path "seqkit_dir/*"
      path "barrnap_rrna_dir/*"
      path checkm_qa
      path gtdb_summary
      path gunc_summary
output:
      path "${params.project}.report.csv"
script:
"""
     bash qc_report_wrapper.sh "${params.project}" seqkit_dir barrnap_rrna_dir $checkm_qa $gtdb_summary $gunc_summary
"""

}

// Detection of chimerism and contamination
 process GUNC {
    tag "${params.project}"

    container params.docker_container_gunc

    publishDir "${params.outdir}/${params.project}/06_GUNC"

input:
      path "genomes/*"

output:
      path "gunc_results/GUNC.progenomes_2.1.maxCSS_level.tsv", emit: gunc_out_ch
script:

"""
  mkdir gunc_results
  gunc run -t 16 --input_dir genomes -r ${params.gunc_db} --file_suffix .${params.ext} -o gunc_results
"""

}

// Generate whole genome
process kSNP4_tree {
    tag "${params.project}"

    errorStrategy 'ignore'
    container params.docker_container_tree
    publishDir "${params.outdir}/${params.project}/07_TREE"

input:
    path "genomes"

output:
    file "kSNP4_output/tree.SNPs_all.ML.tre"

script:
"""
  mkdir kSNP4_output
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

     publishDir "${params.outdir}/${params.project}/07_TREE", mode:'copy', pattern: '*.pdf'

     input:
     file tree

     output:
     file "*.pdf"

     script:
     """
     figtree -graphic PDF "$tree" "${params.project}".bestTree.pdf
     """
 }

