fastq_files_ch = Channel.fromFilePairs(params.fastq_path, checkIfExists: true)

process hisat2_aligment {
	
	publishDir("${params.alignment_out}", mode: 'copy')

	input:
		tuple val(sample_id), path(fastq)

	output:
		path("*.bam")                   , emit: bam
    	path("*.log")                   , emit: summary

	script:

	"""
	hisat2 -q -p $params.threads -x $params.hisat2_index -1 ${fastq[0]} -2 ${fastq[1]} --summary-file ${sample_id}.log | samtools sort -O BAM -o ${sample_id}.bam -

	"""
}


workflow {
	
	hisat2_aligment(fastq_files_ch)
}