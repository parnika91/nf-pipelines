nextflow.enable.dsl=2

import groovy.csv.CsvParser

params.samplesheet = params.samplesheet ?: 'data/samplesheet.csv'

/*
 * Channel: read sample sheet
 */
Channel
  .fromPath(params.samplesheet)
  .map { file(it) }
  .splitCsv(header:true)
  .map { row ->
    // Expect columns: sample, fastq1, fastq2, condition
    def sample = row.sample
    def r1 = file(row.fastq1)
    def r2 = row.fastq2 ? file(row.fastq2) : null
    def cond = row.condition
    [ sample, r1, r2, cond ]
  }
  .set { SAMPLES }

/*
 * STAR alignment (handles PE or SE)
 */
process STAR_ALIGN {
  tag "${sample}"
  publishDir "${params.outdir}/bam", mode:'copy'

  input:
  tuple val(sample), path(r1), path(r2), val(cond)

  output:
  tuple val(sample), path("${sample}.Aligned.sortedByCoord.out.bam"), val(cond)

  when:
  r1

  script:
  def pe = r2 ? true : false
  """
  STAR \\
    --runThreadN ${task.cpus} \\
    --genomeDir ${params.star_index} \\
    --readFilesIn ${r1} ${pe ? r2 : ''} \\
    --readFilesCommand zcat \\
    --outFileNamePrefix ${sample}. \\
    --outSAMtype BAM SortedByCoordinate \\
    --quantMode GeneCounts
  """
}

/*
 * featureCounts on all BAMs together to make one counts matrix
 */
process FEATURECOUNTS {
  tag "featureCounts_all"
  publishDir "${params.outdir}/counts", mode:'copy'

  input:
  tuple val(sample), path(bam), val(cond) from STAR_ALIGN.collect()

  output:
  path "gene_counts.txt"
  path "gene_counts.summary"

  script:
  // list all BAMs and a sample-to-condition table
  def bam_list = bam.collect { it.getName() }.join(' ')
  def sample_names = sample.join('\\n')
  def conds = cond.join('\\n')
  """
  ls -1 *.bam > bam.list

  featureCounts \\
    -T ${task.cpus} \\
    -a ${params.gtf} \\
    -o gene_counts.txt \\
    -g gene_id \\
    -t exon \\
    -s ${params.strandedness} \\
    \$(cat bam.list)

  # build metadata table expected by edgeR script
  paste -d ',' <(echo -e "${sample_names}") <(echo -e "${conds}") > sample_metadata.csv
  """
}

/*
 * edgeR differential expression
 */
process EDGER {
  tag "edgeR_DGE"
  publishDir "${params.outdir}/edgeR", mode:'copy'

  input:
  path "gene_counts.txt" from FEATURECOUNTS.out
  path "sample_metadata.csv" from FEATURECOUNTS.out

  output:
  path "edgeR_DE_results.csv"
  path "edgeR_MAplot.png"
  path "edgeR_MDS.png"

  script:
  """
  Rscript ${projectDir}/bin/edgeR_dge.R \\
    --counts gene_counts.txt \\
    --meta sample_metadata.csv \\
    --outdir .
  """
}

