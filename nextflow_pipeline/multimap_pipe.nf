nextflow.enable.dsl=2

params.data="${projectDir}/../../nextflow_pd/output/fastqc/*R{1,2}*.fastq.gz"

output_dir = "${baseDir}/output"

include { get_packages } from '../modules/get_packages'
include { genome_download } from '../modules/gencode_genome_download'
include { fastp } from '../modules/fastp'
include { Star_genome_gen as star_genome_gen } from '../modules/Star_genome_gen'
include { STAR_multimapping as star_multi_map } from '../modules/STAR_multimapping.nf'
include { sam_sort_index } from '../modules/sam_sort_index'

workflow {
data=Channel.fromFilePairs("${params.data}")
output_dir=Channel.value("${baseDir}/output")
genome_download()
get_packages()
fastp(data, get_packages.out.pack_done_val)
star_genome_gen(genome_download.out.fasta, genome_download.out.gtf)
star_multi_map(fastp.out.reads, star_genome_gen.out.gdir_val)
sam_sort_index(data)
}
