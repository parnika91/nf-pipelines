nextflow.enable.dsl=2

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
  .view()
  .set { SAMPLES }

// Add this validation immediately after building SAMPLES
SAMPLES = SAMPLES.map { s, r1, r2, c ->
    assert r1.exists(), "FASTQ R1 not found for sample ${s}: ${r1}"
    if (r2) assert r2.exists(), "FASTQ R2 not found for sample ${s}: ${r2}"
    [s, r1, r2, c]
}

// Optional: debug print to check what's inside SAMPLES
SAMPLES.view { s, r1, r2, c ->
    "[DEBUG] sample=${s} | r1=${r1} | r2=${r2 ?: '-'} | cond=${c}"
}

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
  echo ">>> [STAR_ALIGN] START sample=${sample} pe=${pe} at \$(date)"
  echo "PWD=\$PWD"
  echo "R1=${r1}"
  echo "R2=${pe ? r2 : '-'}"

  ulimit -n 5000

  STAR \\
    --runThreadN ${task.cpus} \\
    --genomeDir ${params.star_index} \\
    --readFilesIn ${r1} ${pe ? r2 : ''} \\
    --readFilesCommand zcat \\
    --outFileNamePrefix ${sample}. \\
    --outSAMtype BAM SortedByCoordinate \\
    #--limitGenomeGenerateRAM 210000000000

  echo ">>> [STAR_ALIGN] DONE sample=${sample} at \$(date)
  """
}

/*
 * featureCounts on all BAMs together to make one counts matrix
 */
process FEATURECOUNTS {
  tag "featureCounts_all"
  publishDir "${params.outdir}/counts", mode:'copy'

  input:
  tuple val(sample), path(bam), val(cond)

  output:
  path "gene_counts.txt"
  path "sample_metadata.csv"

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
  path "gene_counts.txt"
  path "sample_metadata.csv"

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


workflow {
  // Run STAR per sample
    star_out = STAR_ALIGN(SAMPLES)

    // Collect all BAMs and associated info
    fc_in = star_out.collect()

    // Call FEATURECOUNTS with fc_in
    fc_out = FEATURECOUNTS(fc_in)

    // Pass its outputs to EDGER
    EDGER(fc_out)
}