nextflow.enable.dsl=2

params.samplesheet   = params.samplesheet   ?: 'samplesheet.csv'
params.genome_fasta  = params.genome_fasta  ?: null
params.gtf           = params.gtf           ?: null
params.star_index    = params.star_index    ?: null
params.outdir        = params.outdir        ?: 'results'
params.strandedness  = params.strandedness  ?: 'unstranded'  // unstranded | fr-firststrand | fr-secondstrand

// -------------------------
// Input: read the samplesheet
// -------------------------
Channel
  .fromPath(params.samplesheet)
  .splitCsv(header:true)
  .map { row ->
    def r1 = file(row.fastq1)
    def r2 = file(row.fastq2)
    }
  .set { SAMPLES }


// -------------------------
// STAR alignment per sample
// -------------------------
process STAR_ALIGN {
  tag { id }
  publishDir "${params.outdir}/alignments", mode: 'copy'
  cpus 8
  memory '24 GB'
  input:
    tuple val(id), val(cond), path(r1), path(r2)
    // replicate index path for each sample
    each path(index_dir) from INDEX_DIR
  output:
    tuple val(id), val(cond),
          path("${id}.ReadsPerGene.out.tab"),
          path("${id}.Aligned.sortedByCoord.out.bam") \
      into ALIGN_OUT
  """
  STAR \
    --runThreadN ${task.cpus} \
    --genomeDir ${index_dir} \
    --readFilesIn ${r1} ${r2} \
    --readFilesCommand zcat \
    --outFileNamePrefix ${id}. \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts
  mv ${id}.ReadsPerGene.out.tab ${id}.ReadsPerGene.out.tab
  mv ${id}.Aligned.sortedByCoord.out.bam ${id}.Aligned.sortedByCoord.out.bam
  """
}

workflow {
  take: SAMPLES
  main:
    STAR_ALIGN(SAMPLES, INDEX_DIR)
    //MERGE_COUNTS(ALIGN_OUT)
    //EDGER(MERGE_COUNTS.out)
}
