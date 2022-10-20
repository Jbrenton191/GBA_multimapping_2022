process sam_sort_index {

myDir = file("${projectDir}/output/Samtools")
myDir.mkdir()

publishDir "${projectDir}/output/Samtools", mode: 'copy', overwrite: true

input:
path(data)

output:
path("*.bai"), emit: bam_indexes
path("*.bam"), emit: sorted_bams
val(bam_dir), emit: bam_dir

script:
bam_dir="${projectDir}/output/Samtools"
"""
name=`echo "$data" | sed 's/_mapped.*//g'`
samtools sort -m 1000000000 $data -o ./\${name}_Aligned.sortedBysamtools.out.bam

samtools index ./\${name}_Aligned.sortedBysamtools.out.bam
"""
}
