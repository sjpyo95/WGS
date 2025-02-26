#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.inputdir
params.outputdir
params.ref_genome
params.window_size
params.trim
params.meandepth
params.vcf_dir
params.adaptor_seq
params.adaptor_seq2

println "Input dir : ${params.inputdir}"
println "Reference : ${params.ref_genome}"
println "Trimming  : ${params.trim}"
println "WindowSize: ${params.window_size}"

Channel
	.fromFilePairs("${params.inputdir}/*_1.fastq.gz", flat: true, size: 2)
	.map{ tuple, files ->
		def sample_name = tuple.getKey().name.split('_')[0]
		tuple.setKey(sample_name)
		return tuple
	}
	.set { reads_ch }

known_sites_ch = Channel
	.fromPath("${params.vcf_dir}/*.vcf.gz")
	.toList()
	.map { vcfs -> vcfs.collect{ "--known-sites " + it }.join(" ") }

process FASTQC {
	tag "${sample}"
	input:
	tuple val(sample), path(reads)
	output:
	path("${sample}_1_fastqc.zip"), path("${sample}_2_fastqc.zip")

	"""
	fastqc ${reads} -o . -t 14
	"""
}

process TRIM_GALORE {
	when:
	params.trim == true

	tag "${sample}"
	input:
	tuple val(sample), path(reads)
	output:
	tuple val(sample), path("${sample}_1_val_1.fq.gz"), path("${sample}_2_val_2.fq.gz")

	"""
	trim_galore --paired --adapter ${params.adaptor_seq} --adapter2 ${params.adaptor_seq2} -q 20 --cores 14 ${reads[0]} ${reads[1]}
}
	"""

reads_for_align = params.trim \
	? TRIM_GALORE.out.map { sample, r1, r2 -> tuple(sample, [r1, r2]) }
	: reads_ch.map { sample, r -> tuple(sample,r) }

process FASTQC_TRIM {
	when:
	params.trim == true

	tag "${sample}"
	input:
	tuple val(sample), path(reads)
	output:
	path("${sample}_1_val_1.fastqc.zip"), path("${sample}_2_val_2.fastqc.zip")

	"""
	fastqc ${reads} -o . -t 14
	"""
}

process BWA_MEM {
	tag "${sample}")
	input:
	tuple val(sample), path(reads) from reads_for_align
	output:
	path("${sample}.sorted.bam")

	"""
	bwa mem -t 14 ${params.ref_genome} ${reads[0]} ${reads[1]} | \
	samtools view -@ 14 -Sb - | \
	samtools sort -o ${sample}.sorted.bam -
	"""
}

process INDEX_BAM {
	tag "${sample}"
	input:
	tuple val(sample), path(bam)
	output:
	path("${bam}.bai")

	"""
	samtools index ${bam}
	"""
}

process COVERAGE {
	tag "${sample}"
	input:
	tuple(sample), path(bam) from BWA_MEM.out
	path(bai) from INDEX_BAM.out
	output:
	path("${bam.baseName}._depth.bed")

	"""
	samtools depth -a ${bam} > ${bam.baseName}_depth.bed
	"""
}

process CAL_DEPTH {
	tag "${sample}"
	input:
	tuple(sample), path(bam) from BWA_MEM.out
	path(bai) from INDEX_BAM.out
	output:
	path("${bam.baseName}_windowed_depth.${(params.window_size/1000 as int)}kbp.bed")

	script:
	"""
	python scripts/window_depth.py ${bam} ${params.window_size} \
		${bam.baseName}_windowed_depth.${(params.window_size/1000 as int)}kbp.bed
	"""
}

process PLOT_DEPTH {
	input:
	path(depth_files) from CAL_DEPTH.out.collect()

	output:
	path("${(params.window_size/1000 as int)}kbp.depth_plot.png"),
	path("${(params.window_size/1000 as int)}kbp.zscore_plot.png")

	script:
	"""
	Rscript scripts/depth_plot.R . \
		${(params.window_size/1000 as int)}kbp.depth_plot.png \
		${(params.window_size/1000 as int)}kbp.zscore_plot.png \
		${params.window_size} ${params.meandepth}
	"""
}

process MARK_DUPLICATES {
	tag "${sample}"
	input:
	tuple val(sample), path(bam) from BWA_MEM.out
	path(bai) from INDEX_BAM.out

	output:
	tuple val(sample),
	path("${bam.baseName}.dedup.bam"),
	path("${bam.baseName}.dedup.bai"),
	path("${bam.baseName}.dedup.metrics.txt")

	"""
	gatk MarkDuplicates \
		-I ${bam} \
		-O ${bam.baseName}.dedup.bam \
		-M ${bam.baseName}.dedup.metrics.txt \
		--CREATE_INDEX true
	"""
}

samples_ch = reads_ch.map { s, _ -> s }.unique()

fastqc_zip_ch = FASTQC.out
	.groupTuple()
	.map { files ->
		def sname = files[0].name.split('_')[0]
		return jtuple(sname, files)
	}

dedup_bam_ch = MARK_DUPLICATES.out
	.groupTuple()
	.map { group ->
		def sname = group[0][0] //first sample tuple
		def bamFile = group.find { it[1].name.endsWith('.dedup.bam') }[1]
		def baiFile = group.find { it[1].name.endWith('.dedup.bai') }[1]
		def metFile = group.find { it[1].name.endsWith('.dedup.metrics.txt') }[1]
		return tuple(sname, bamFile, baiFile, metFile)
	}

process QC_STAT {
	tag "${sample}"
	input:
	tuple val(sample), val(fqzips), val(bam), val(bai), val(met)
	output:
	path("${sample}_stat.txt")

	script:
	"""
	python scripts/qc_stat.py \
		${fqzips[0]} ${fqzips[1]} \
		${bam} \
		${params.ref_genome}.fai \
		${sample}_stat.txt
	"""
}

samples_ch
	.join(fastqc_zip_ch)
	.join(dedup_bam_ch)
	.map { s, fc, dedup ->
		def sampleName = s
		def fqzips = fc[1]
		def bamFile = dedup[1]
		def baiFile = dedup[2]
		def metFile = dedup[3]
		return tuple(sampleName, fqzips, bamFile, baiFile, metFile)
	}
	.set { qc_input_ch }

QC_STAT.inChannel(qc_input_ch)

process AGGREGATE_STAT {
	input:
	path(stat_files) from QC_STAT.out.collect()
	output:
	path("qc_stat.txt")

	script:
	"""
	python scripts/aggregate_stat.py . qc_stat.txt
	"""
}


process HAPLOTYPECALLER {
	tag "${sample}"

	input:
	tuple val(sample), path(bam), path(bai), path(met) from dedup_bam_ch
	output:
	path("${sample}.g.vcf.gz")

	"""
	gatk HaplotypeCaller \
		-R ${params.ref_genome} \
		-I ${bam} \
		-O ${sample}.g.vcf.gz \
		-ERC GVCF
	"""
}

process GENOTYPEGVCFS {
	input:
	path vcf_files
	output:
	path("combined_g.vcf.gz")

	script:
	"""
	gatk GenotypeGVCFs -R ${params.ref_genome} \
		${vcf_files.collect{ "--variant " + it }.join(" ")} \
		-O combined_g.vcf.gz
	"""
}

workflow {
	out_fastqc = FASTQC(reads_ch)

	if (params.trim) {
		out_trim = TRIM_GALORE(reads_ch)
		FASTQC_TRIM(out_trim)
		reads_for_align = out_trim.map { sample, r1, r2 -> tuple(sample, [r1, r2]) }
	} else {
		reads_for_align = reads_ch
	}

	out_bwa = BWA_MEM(reads_for_align)
	

