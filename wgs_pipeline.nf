params.fastq_dir = '/path/to/fastq_files/'
params.outputdir = '/path/to/output/directory/'
params.ref_genome = '/path/to/ref.genome.fa'
params.threads = 4
params.min_threads = 4
params.start_step = 'fastqc'
params.end_step = 'calculate_qc_metrics'

// Helper function to check step range
def step_in_range(step) {
	def steps = ['fastqc', 'trimgalore', 'bwa_alignment', 'create_sample_list', 'calculate_qc_metrics']
	def start_idx = steps.indexOf(params.start_step)
	def end_idx = steps.indexOf(params.end_step)
	def current_idx = steps.indexOf(step)
	return current_idx >= start_idx && current_idx <= end_idx
}

// Function to calculate threads per sample
def calculate_threads_per_sample(sample_count) {
	def threads_per_sample = (params.threads / sample_count).intValue()
	return Math.max(threads_per_sample, params.min_threads)
}

// Function to calculate parallel samples to run
def calculate_parallel_samples(sample_count) {
	def threads_per_sample = calculate_threads_per_sample(sample_count)
	def parallel_samples = (params.threads / threads_per_sample).intValue()
	return Math.min(parallel_samples, sample_count) // Cannot run more than the sample count
}


// Flags to confrol which steps to run
params.run_fastqc = true
params.run_trimgalore = false
params.run_alignment = true
params.run_sample_list = true
params.run_qc_metrics = true

process fastqc {
	tag "${sample_id}"

	input:
	tuple val(sample_id), path(fastq1), path(fastq2)

	output:
	path "${params.output_dir}fastqc/${sample_id}_fastqc.zip"

	when:
	step_in_range('fastqc')

	script:
	"""
	python fastQC.py --fastqs ${fastq1} ${fastq2} --outputdir ${params.outputdir}/fastqc --threads ${task.cpus}
	"""
}

process trimgalore {
	tag "${sample_id}"

	input:
	tuple val(sample_id), path(fastq1), path(fastq2)

	output:
	path "${params.outputdir}/trimmed_fastq/${sample_id}_R1_trimmed.fastq.gz", path "${params.outputdir}/trimmed_fastq/${sample_id}_R2_trimmed.fastq.gz"

	when:
	params.run_trimgalore && step_in_range('trimgalore')

	script:
	"""
	python trimgalore.py --fastq1 ${fastq1} --fastq2 ${fastq2} --outputdir ${params.outputdir}/trimmed_fastq/ --threads ${task.cpus}
	"""
}

process bwa_alignment {
	tag "${sample_id}"

	input:
	tuple val(sample_id), path(fastq1), path(fastq2)

	output:
	path "${params.outputdir}/alignment/${sample_id}.sorted.bam"

	when:
	step_in_range('bwa_alignment')

	script:
	"""
	python bwa_alignment.py --fastq1 ${fastq1} --fastq2 ${fastq2} --ref_genome ${params.ref_genome} --outputdir ${params.outputdir}/alignment/ --threads ${task.cpus}
	"""
}

process generate_bed {
	tag "${sample_id}"

	input:
	tuple val(sample_id), path(bam_file)

	output:
	path "${params.output_dir}/beds/${sample_id}.bed"

	script:
	"""
	python generate_bed_from_bam.sh ${bam_file} ${params.output_dir}/beds/
	"""
}
process create_sample_list {
	input:
	path "${params.outputdir}/alignment/*.sorted.bam"

	output:
	path "${params.outputdir}/sample_list.csv"

	when:
	step_in_range('create_sample_list')

	script:
	"""
	ptyhon create_sample_list.py --fastqc_dir ${params.outputdir}/fastqc/ --bam_dir ${params.outputdir}/alignment/ --output_file ${params.outputdir}/sample_list.csv
	"""
}

process calculate_qc_metrics {
	input:
	path "${params.outputdir}/sample_list.csv"

	output:
	path "${params.outputdir}/multi_sample_qc_metrics.csv"

	when:
	step_in_range('calculate_qc_metrics')

	script:
	"""
	python calculate_qc_metrics.py --sample_list ${params.outputdir}/sample_list.csv --ref_genome ${params.ref_genome} --outputdir ${params.outputdir}
	"""
}

workflow {
	// Get all FASTQ files and count the number of samples
	fastq_pairs = Cannel
		.fromFilePairs("${params.fastq_dir}*_R{1,2}.fatq.gz", flat:true)
		.map { sample_id, fastqs -> tuple(sample_id, fastqs[0], fastqs[1]) }

	sample_count = fastq_pairs.count()
	
	// Calculate threads per sample and parallel samples to run
	def threads_per_sample = calculate_threads_per_sample(sample_count)
	def parallel_samples = calculate_parallel_samples(sample_count)

	fastq_pars
		.set { all_fastq_pairs }


	// Step 1: FastQC
	all_fastq_pairs
		.buffer(parallel_samples)
		| fastqc

	// Step 2: Trim Galore
	if (params.run_trimgalore && step_in_range('trimgalore')) {
		all_fastq_pairs
			.buffer(parallel_samples)
			| trimgalore
	}

	// Step 3: BWA Alignment
	all_fastq_pairs
		.buffer(parallel_samples)
		| bwa_alignment
	
	// Step 4: BAM to BED conversion (Depth calculation)
	if (step_in_range('generate_bed')) {
		aligned_bams
			.set { all_bams }
		all_bams
			| generate_bed
	}

	// Step 5: Create Sample List and Calculate QC Metrics (executed once after alignment for all samples)
	if (step_in_range('create_sample_list')) {
		bwa_alignment.out.collect()
		| create_sample_list
	}

	if (step_in_range('calculate_qc_metrics')) {
		create_sample_list.out.collect()
		| calculate_qc_metrics
	}
}
