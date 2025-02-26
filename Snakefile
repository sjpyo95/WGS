#!~/anaconda3/envs/ctcells/bin/python
# WGS Pipeline
# FastQC to Germline Variant Calling Steps
# Config Path: ./config/config.yaml

import glob
import re
import yaml
import os

configfile: 'config/config.yaml'


inputdir = config['inputdir']
outputdir = config['outputdir']
ref_genome = config['ref_genome']
window_size = config['window_size']
window_size_kbp = str(int(int(window_size)/1000))
trimming = config['trim']

fastqs= glob.glob(f"{inputdir}/*.fastq.gz")
print(fastqs)
sample_to_reads = {}
meandepth = config['meandepth']

for fastq in fastqs:
	filename = os.path.basename(fastq)
	sample_name = filename.split("_")[0]
	read = filename.split("_")[1][0]
	
	if sample_name not in sample_to_reads:
		sample_to_reads[sample_name] = {}
	sample_to_reads[sample_name][read] = fastq

SAMPLES = list(sample_to_reads.keys())

def getFastqfiles(wildcards):
	if trimming:
		return [
			f"{outputdir}/trim/{{sample}}_1_val_1.fq.gz",
			f"{outputdir}/trim/{{sample}}_2_val_2.fq.gz"
		]
	else:
		return [
			f"{inputdir}/{{sample}}_1.fastq.gz",
			f"{inputdir}/{{sample}}_2.fastq.gz"
		]

def getFastqczip(wildcards):
	if trimming:
		return [
			f"{outputdir}/fastqc_trim/{wildcards.sample}_1_val_1_fastqc.zip",
			f"{outputdir}/fastqc_trim/{wildcards.sample}_2_val_2_fastqc.zip"
		]
	else:
		return [
			f"{outputdir}/fastqc/{wildcards.sample}_1_fastqc.zip",
			f"{outputdir}/fastqc/{wildcards.sample}_2_fastqc.zip"
		]

def get_known_sites_params(wildcards):
	vcfs = filter(lambda x: '.vcf' in x and not '.tbi' in x, os.listdir(config['vcf_dir']))
	known_sites = ['--known-sites ' + config['vcf_dir'] + '/' + x for x in vcfs]
	return ' '.join(known_sites)

# ------------ Final results ------------ #
rule all:
	input:
		expand(f"{outputdir}/fastqc/{{sample}}_1_fastqc.zip", sample=SAMPLES),
		expand(f"{outputdir}/fastqc_trim/{{sample}}_1_val_1_fastqc.zip" if trimming else f"{outputdir}/fastqc/{{sample}}_1_fastqc.zip", sample=SAMPLES),
		expand(f"{outputdir}/fastqc_trim/{{sample}}_2_val_2_fastqc.zip" if trimming else f"{outputdir}/fastqc/{{sample}}_2_fastqc.zip", sample=SAMPLES),
		expand(f"{outputdir}/depth/{window_size_kbp}kbp/{{sample}}_windowed_depth.{window_size_kbp}kbp.bed", sample=SAMPLES),
		f"{outputdir}/depth/{window_size_kbp}kbp/{window_size_kbp}kbp.zscore_plot.png",
		f"{outputdir}/qc/qc_stat.txt",
		expand(f"{outputdir}/variants/{{sample}}.g.vcf.gz", sample=SAMPLES)


# ------------ FastQC ------------ #
rule fastqc:
	input:
		R1=f"{inputdir}/{{sample}}_1.fastq.gz",
		R2=f"{inputdir}/{{sample}}_2.fastq.gz"
	output:
		zip1=f"{outputdir}/fastqc/{{sample}}_1_fastqc.zip",
		zip2=f"{outputdir}/fastqc/{{sample}}_2_fastqc.zip"
	threads:
		14
	log:
		f'{outputdir}/fastqc/logs/{{sample}}_fastqc.log'
	shell:
		'''
		fastqc {input.R1} {input.R2} -o {outputdir}/fastqc -t {threads} >> {log} 2>&1
		'''
# ------------ Trim Galore (optional) ------------ #
rule trim:
	input:
		R1=f"{inputdir}/{{sample}}_1.fastq.gz",
		R2=f"{inputdir}/{{sample}}_2.fastq.gz"
	output:
		tR1=f"{outputdir}/trim/{{sample}}_1_val_1.fq.gz",
		tR2=f"{outputdir}/trim/{{sample}}_2_val_2.fq.gz"
	params:
		adp1=config['adaptor_seq'],
		adp2=config['adaptor_seq2'],
		q='20'
	log:
		f"{outputdir}/trim/logs/{{sample}}_trim.log"
	threads:
		14
	shell:
		'''
		trim_galore --paired --adapter {params.adp1} --adapter2 {params.adp2} --o {outputdir}/trim/ -q {params.q} --cores {threads} {input.R1} {input.R2} >> {log} 2>&1
		'''

rule fastqc_trim:
	input:
		tR1=f"{outputdir}/trim/{{sample}}_1_val_1.fq.gz",
		tR2=f"{outputdir}/trim/{{sample}}_2_val_2.fq.gz"
	
	output:
		tzip1=f"{outputdir}/fastqc_trim/{{sample}}_1_val_1_fastqc.zip",
		tzip2=f"{outputdir}/fastqc_trim/{{sample}}_2_val_2_fastqc.zip",
	threads:
		14
	log:
		f"{outputdir}/fastqc_trim/logs/{{sample}}_fastqc_trim.log"
	shell:
		'''
		fastqc {input.tR1} {input.tR2} -o {outputdir}/fastqc_trim -t {threads} >> {log} 2>&1
		'''

		
# ------------ BWA-MEM Alignment ------------ #
rule bwa:
	input:
		fastqs=getFastqfiles,
		ref=ref_genome
	output:
		temp(f"{outputdir}/bwa/{{sample}}.bam")
	threads:
		14
	shell:
		'''
		bwa mem -t {threads} {input.ref} {input.fastqs[0]} {input.fastqs[1]} | samtools view -@ {threads} -Sb - > {output}
		'''

# ------------ Sort BAM ------------ #
rule sort_bam:
	input:
		f"{outputdir}/bwa/{{sample}}.bam"
	output:
		f"{outputdir}/bwa/{{sample}}.sorted.bam"
	shell:
		'''
		samtools sort -o {output} {input}
		'''

# ------------ Index Sorted BAM File ------------ #
rule index_bam:
    input:
        f"{outputdir}/bwa/{{sample}}.sorted.bam"
    output:
        f"{outputdir}/bwa/{{sample}}.sorted.bam.bai"
    shell:
        "samtools index {input}"

# 
rule coverage:
	input:
		bam=f"{outputdir}/bwa/{{sample}}.sorted.bam",
		bai=f"{outputdir}/bwa/{{sample}}.sorted.bam.bai"
	output:
		f"{outputdir}/depth/{{sample}}_depth.bed"
	shell:
		'''
		samtools depth -a {input.bam} > {output}
		'''


# ------------ Window Depth BED Files ------------ #
rule cal_depth:
	input:
		bam=f"{outputdir}/bwa/{{sample}}.sorted.bam",
		bai=f"{outputdir}/bwa/{{sample}}.sorted.bam.bai"
	output:
		depth=f"{outputdir}/depth/{window_size_kbp}kbp/{{sample}}_windowed_depth.{window_size_kbp}kbp.bed"
	threads: 
		14
	shell:
		'''
		python scripts/window_depth.py {input.bam} {window_size} {output.depth}
		'''
# ------------ Genome Wide Depth Coverage Polt ------------ #
rule plot_depth:
	input:
		depth=expand(f"{outputdir}/depth/{window_size_kbp}kbp/{{sample}}_windowed_depth.{window_size_kbp}kbp.bed", sample=SAMPLES)
	output:
		depth_plot=f"{outputdir}/depth/{window_size_kbp}kbp/{window_size_kbp}kbp.depth_plot.png",
		zscore_plot=f"{outputdir}/depth/{window_size_kbp}kbp/{window_size_kbp}kbp.zscore_plot.png"
	params:
		inputdir=f"{outputdir}/depth/{window_size_kbp}kbp"
	log:
		f"{outputdir}/depth/logs/{window_size_kbp}kbp_plot_depth.log"
	shell:
		'''
		Rscript scripts/depth_plot.R {params.inputdir} {output.depth_plot} {output.zscore_plot} {window_size} {meandepth} >> {log} 2>&1
		'''

# ------------ GATK4 Mark Duplicates ------------ #
rule mark_duplicates:
	input:
		bam=f"{outputdir}/bwa/{{sample}}.sorted.bam",
		bai=f"{outputdir}/bwa/{{sample}}.sorted.bam.bai"
	output:
		bam=f"{outputdir}/gatk/{{sample}}.dedup.bam",
		bai=f"{outputdir}/gatk/{{sample}}.dedup.bai",
		metrics=f"{outputdir}/gatk/{{sample}}.dedup.metrics.txt"
	log:
		f"{outputdir}/gatk/logs/{{sample}}_markDuplicates.log"
	shell:
		'''
		gatk MarkDuplicates \
			-I {input.bam} \
			-O {output.bam} \
			-M {output.metrics} \
			--CREATE_INDEX true >> {log} 2>&1
		'''
# ------------ make QC Stat metrics ------------- #
			
rule qc_stat:
	input:
		zips=getFastqczip,
		bam=f"{outputdir}/gatk/{{sample}}.dedup.bam",
		bai=f"{outputdir}/gatk/{{sample}}.dedup.bai",
		metrics=f"{outputdir}/gatk/{{sample}}.dedup.metrics.txt"
	output:
		f"{outputdir}/qc/{{sample}}_stat.txt"
	params:
		fai=ref_genome+'.fai'
	shell:
		'''
		python scripts/qc_stat.py {input.zips[0]} {input.zips[1]} {input.bam} {params.fai} {output}
		'''
# ------------ Merge stat metrics ------------ #
rule aggregate_stat:
	input:
		expand(f"{outputdir}/qc/{{sample}}_stat.txt", sample=SAMPLES)
	output:
		f'{outputdir}/qc/qc_stat.txt'
	run:
		headers = [
			'Sample ID', 
			'% > Fastq Q20', 
			'% > Fastq Q30',
			'Total reads',
			'Read length (bp)',
			'Total yield (Mbp)',
			'Reference size (Mbp)',
			'Throughput mean depth (X)',
			'De-duplicated reads',
			'De-duplicated reads % (out of total reads)',
			'Mappable reads (reads mapped to human genome)',
			'Mappable reads % (out of de-duplicated reads)',
			'Mappable yield (Mbp)',
			'Mappable mean depth (X)',
			'Mean Depth (dedup) (X)',
			'% >= 1X coverage',
			'% >= 5X coverage',
			'% >= 10X coverage',
			'% >= 15X coverage',
			'% >= 20X coverage',
			'% >= 30X coverage',
		]
		
		with open(output[0], 'w') as outfile:
			outfile.write('\t'.join(headers))

			for samplefile in input:
				with open(samplefile, 'r') as infile:
					for line in infile:
						outfile.write(line)

# ------------ Add Read Groups ------------ #
rule add_readgroups:
	input:
		bam=f"{outputdir}/gatk/{{sample}}.dedup.bam",
		bai=f"{outputdir}/gatk/{{sample}}.dedup.bai"
	output:
		bam=f"{outputdir}/gatk/{{sample}}.rg.bam"
	params:
		rgid='1',
		rglb='lib1',
		rgpl='illumina',
		rgpu='ctcells',
		rgsm=lambda wildcards: wildcards.sample
	shell:
		'''
		gatk AddOrReplaceReadGroups \
			-I {input.bam} \
			-O {output.bam} \
			-RGID {params.rgid} \
			-RGLB {params.rglb} \
			-RGPL {params.rgpl} \
			-RGPU {params.rgpu} \
			-RGSM {params.rgsm}
		'''

rule rg_index:
	input:
		bam=f"{outputdir}/gatk/{{sample}}.rg.bam"
	output:
		bai=f"{outputdir}/gatk/{{sample}}.rg.bam.bai"
	shell:
		'''
		samtools index {input.bam}
		'''

# ------------ GATK4 BaseRecalibroator ------------ #
rule base_recal:
	input:
		bam=f"{outputdir}/gatk/{{sample}}.rg.bam",
		bai=f"{outputdir}/gatk/{{sample}}.rg.bam.bai"
	output:
		recal_table=f"{outputdir}/gatk/{{sample}}.recal_data.table"
	params:
		ref=ref_genome,
		known_sites=get_known_sites_params
	shell:
		'''
		gatk BaseRecalibrator \
			-I {input.bam} \
			-R {params.ref} \
			{params.known_sites} \
			-O {output.recal_table}
		'''

# ------------ GATK4 ApplyBQSR ------------ #
rule apply_bqsr:
	input:
		recal_table=f"{outputdir}/gatk/{{sample}}.recal_data.table",
		ref=ref_genome
	output:
		bam=f"{outputdir}/gatk/{{sample}}.recal.bam",
		bai=f"{outputdir}/gatk/{{sample}}.recal.bam.bai"
	params:
		bam=f"{outputdir}/gatk/{{sample}}.rg.bam",
		bai=f"{outputdir}/gatk/{{sample}}.rg.bai",

	shell:
		'''
		gatk ApplyBQSR \
			-R {input.ref} \
			-I {params.bam} \
			--bqsr-recal-file {input.recal_table} \
			-O {output.bam}
		'''

# ------------ GATK4 HaplotypeCaller ------------ #
rule haplotypecaller:
	input:
		bam=f"{outputdir}/gatk/{{sample}}.recal.bam",
		bai=f"{outputdir}/gatk/{{sample}}.recal.bam.bai",
		ref=ref_genome
	output:
		vcf=f"{outputdir}/variants/{{sample}}.g.vcf.gz"
	shell:
		'''
		gatk HaplotypeCaller \
			-R {input.ref} \
			-I {input.bam} \
			-O {output.vcf} \
			-ERC GVCF
		'''

rule genotypeGVCFs:
	input:
		expand(f"{outputdir}/variants/{{sample}}.g.vcf.gz", sample=SAMPLES)
	output:
		vcf = f"{outputdir}/variants/combined.g.vcf.gz"
	params:
		ref=ref_genome
	shell:
		"""
		gatk GenotypeGVCFs \
			-R {params.ref} \
			{(" ").join(["--variant " + v for v in input])} \
			-O {output.vcf}
		"""
