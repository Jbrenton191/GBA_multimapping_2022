process get_packages {
	
	publishDir "${projectDir}", mode: 'copy', overwrite: true

	output:
//	path("*")
	val(fin_val), emit: pack_done_val
	
	script:
	fin_val="Package download didn't fail"
	"""
	Rscript $projectDir/R_scripts/Rpackage_download.R
	"""

// add in (also move yml file into folder): conda env create -f nextflow_env.yml -n nf_env
}
