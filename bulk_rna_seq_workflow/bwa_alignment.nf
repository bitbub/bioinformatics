fastq_files_ch = Channel.fromFilePairs(params.fastq_path, checkIfExists: true)

process bwa_aligment {
	
	publishDir("${params.alignment_out}", mode: 'copy')

	input:
		tuple val(sample_id), path(fastq)

	output:
		path "*"

	script:

	"""
	
	bwa mem -t $params.threads $params.bwa_index ${fastq} | samtools sort -@ $params.threads -O BAM -o ${sample_id}.bam -
	
	"""
}


workflow {
	
	bwa_aligment(fastq_files_ch)
}