process STAR_multimapping {

myDir2 = file("${projectDir}/output/STAR/align")
myDir2.mkdirs()

publishDir "${projectDir}/output/STAR/align", mode: 'copy', overwrite: true

  input:
  tuple val(sampleID), path(reads)
  val(genome_dir)

  output:
	path("*SJ.out.tab"), emit: sj_tabs2
  	val(sj_loc), emit: sj_loc
        path('*.bam'), emit: bams

        tuple val(sampleID), path('*toTranscriptome.out.bam'), optional:true, emit: bam_transcript
        tuple val(sampleID), path('*Aligned.unsort.out.bam'), optional:true, emit: bam_unsorted

        tuple val(sampleID), path('*Log.final.out'), emit: log_final
        tuple val(sampleID), path('*Log.out'), emit: log_out
        tuple val(sampleID), path('*Log.progress.out'), emit: log_progress

  script:
  sj_loc="${projectDir}/output/STAR/align/pre_merge"
  """
  echo ${sampleID}
  echo ${reads[0]}
  echo ${reads[1]}

  STAR --runThreadN 7 \
--genomeDir $genome_dir \
--readFilesIn ${reads[0]}, ${reads[1]} \
--readFilesCommand zcat \
--outFileNamePrefix ${sampleID}_mapped.BAM_ \
--outSAMtype BAM Unsorted \
--outFilterType BySJout \
--outFilterMultimapNmax 10 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverReadLmax 0.04 \
--alignIntronMin 25 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 3 \
--outSAMmultNmax -1
"""
// --outSAMunmapped Within \
}
