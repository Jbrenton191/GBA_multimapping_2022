process Star_genome_gen {

myDir = file("${projectDir}/output/STAR")
myDir.mkdirs()

myDir2 = file("${projectDir}/output/STAR/genome_dir")
myDir2.mkdirs()

publishDir "${projectDir}/output/STAR/genome_dir", mode: 'copy', overwrite: true

    input:
    path(fasta)
    path(gtf)

    output:
    val("${gdir_val}"), emit: gdir_val

    script:
    gdir_val=file("${projectDir}/output/STAR/genome_dir")
//    ref_dir="${baseDir}/output/reference_downloads"
    """
STAR --runThreadN 7 \
--runMode genomeGenerate \
--genomeDir $gdir_val \
--genomeFastaFiles $fasta \
--sjdbGTFfile $gtf \
--sjdbOverhang 99
  """
}
