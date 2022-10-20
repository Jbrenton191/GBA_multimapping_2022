process genome_download {


myDir3 = file("${baseDir}/output/reference_downloads")
myDir3.mkdir()

publishDir "${baseDir}/output/reference_downloads", mode: 'copy', overwrite: true

	output:
	path("*primary_assembly*.fa"), emit: fasta
	path("*.gtf"), emit: gtf
	path("*.transcripts.fa"), emit: transcripts

	script:
	"""
	wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
	wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz

	gunzip GRCh38.primary_assembly.genome.fa.gz
	gunzip gencode.v38.annotation.gtf.gz

	wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.transcripts.fa.gz

	gunzip gencode.v38.transcripts.fa.gz
	"""
}
