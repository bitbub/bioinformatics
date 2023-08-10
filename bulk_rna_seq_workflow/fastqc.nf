#! /usr/bin/env nextflow

fastq_files_ch = Channel.fromPath(params.fastq_path, checkIfExists: true)


process fastqc {
	
	publishDir( "${params.fastqc_out}", mode: 'copy' )

	input:
		path fastq_file

	output:
		path "*"

	script:

		"""
		fastqc -q -t $params.threads ${fastq_file}
		"""
}


workflow {
	fastqc(fastq_files_ch)
}