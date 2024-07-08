#conda activate CHIP
import subprocess
import os, sys
from glob import glob
from util import run_command, create_multi_dir, run_parallel

def variant_calling(sample_id, input_dir, output_dir, tps, reference_genome, known_sites_vcfs):
	'''Process a single sample for variant calling.'''
	known_sites_list = known_sites_vcfs.split(',')
	
	# Find corresponding FASTQ files
	r1_files = glob(os.path.join(input_dir, sample_id, "*_1.fastq.gz"))
	r2_files = glob(os.path.join(input_dir, sample_id, "*_2.fastq.gz"))
	if not r1_files or not r2_files:
		print(f"No FASTQ files found for sample {sample_id}")
		return None
	
	r1 = r1_files[0]
	r2 = r2_files[0]

	# Step 1: Quality Control
	run_command(f"fastqc -t {tps} {r1} {r2} -o {output_dir}/qc")

	# Step 2: Trimming using Trim Galores
	run_command(f"trim_galore --paired --cores {tps} --output_dir {output_dir}/trimmed {r1} {r2}")

	# Step 3: Alignment
	r1_paired = os.path.join(output_dir, 'trimmed', os.path.basename(r1).replace('.fastq.gz', '_val_1.fq.gz'))
	r2_paired = os.path.join(output_dir, 'trimmed', os.path.basename(r2).replace('.fastq.gz', '_val_2.fq.gz'))
	run_command(f"bwa mem -t {tps} {reference_genome} {r1_paired} {r2_paired} > {output_dir}/aligned/{sample_id}.sam")

	# Step 4: Convert SAM to BAM, sort, and index
	run_command(f"samtools view -@ {tps} -bS {output_dir}/aligned/{sample_id}.sam > {output_dir}/aligned/{sample_id}.bam")
	run_command(f"samtools sort -@ {tps} {output_dir}/aligned/{sample_id}.bam -o {output_dir}/sorted/{sample_id}_sorted.bam")
	run_command(f"samtools index {output_dir}/sorted/{sample_id}_sorted.bam")

   # Step 5: Add Read Groups
	run_command(f"picard AddOrReplaceReadGroups I={output_dir}/sorted/{sample_id}_sorted.bam O={output_dir}/sorted/{sample_id}_rg.bam RGID={sample_id} RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM={sample_id}")
	run_command(f"samtools index {output_dir}/sorted/{sample_id}_rg.bam")

	#Step 6: Mark Duplicates
	run_command(f"picard MarkDuplicates I={output_dir}/sorted/{sample_id}_rg.bam O={output_dir}/deduped/{sample_id}_deduped.bam M={output_dir}/deduped/{sample_id}_deduped.metrics CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT TMP_DIR={output_dir}/deduped")

    # Step 7: Base Quality Score Recalibration (BQSR)
	known_sites_args = ' '.join([f"--known-sites {vcf}" for vcf in known_sites_list])
	run_command(f"gatk BaseRecalibrator -I {output_dir}/deduped/{sample_id}_deduped.bam -R {reference_genome} {known_sites_args} -O {output_dir}/bqsr/{sample_id}_recal_data.table")
	run_command(f"gatk ApplyBQSR -R {reference_genome} -I {output_dir}/deduped/{sample_id}_deduped.bam --bqsr-recal-file {output_dir}/bqsr/{sample_id}_recal_data.table -O {output_dir}/bqsr/{sample_id}_recalibrated.bam")
	run_command(f"samtools index {output_dir}/bqsr/{sample_id}_recalibrated.bam")
	
	# Step 8: Variant Calling with HaplotypeCaller (GVCF mode)
   gvcf_path = f"{output_dir}/variants/{sample_id}.g.vcf.gz"
	run_command(f"gatk HaplotypeCaller -R {reference_genome} -I {output_dir}/bqsr/{sample_id}_recalibrated.bam -O {gvcf_path} -ERC GVCF --native-pair-hmm-threads {tps}")
	
	return gvcf_path

def main(input_dir, output_dir, reference_genome, known_sites_vcfs, total_threads, samples_per_run):
	# Split the comma-separated known sites VCFs into a list
	known_sites_list = known_sites_vcfs.split(',')

	# Create output directories
	directories = ["qc", "trimmed", "aligned", "sorted", "deduped", "bqsr", "variants", "final_variants", "annotated"]
	create_multi_dir(output_dir, directories)

	# Process each sample directory in the input directory
	sample_dirs = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]
	sample_dirs.sort()
	sample_dirs = sample_dirs[0:3]
	
	# Run the pipeline in parallel
	gvcf_list = run_parallel(sample_dirs, input_dir, output_dir, total_threads, samples_per_run, variant_calling, reference_genome, known_sites_vcfs)
		    
	# Step 9: Joint Genotyping
	# Create a text file with the list of GVCF files
	gvcf_list_path = os.path.join(output_dir, "variants/gvcf_list.txt")
	with open(gvcf_list_path, 'w') as f:
		for gvcf in gvcf_list:
			f.write(gvcf+'\n')
	
	# Combine GVCFs
	combined_gvcf = os.path.join(output_dir, 'variants/combined.g.vcf.gz')
	gvcf_args = ' '.join([f"--variant {gvcf}" for gvcf in gvcf_list])
	run_command(f"gatk CombineGVCFs -R {reference_genome} {gvcf_args} -O {combined_gvcf}")

	# Genotype GVCFs
	raw_vcf = os.path.join(output_dir, "variants/raw_variants.vcf.gz")
	run_command(f"gatk GenotypeGVCFs -R {reference_genome} -V {combined_gvcf} -O {raw_vcf}")

	# Step 9: Varitant Filtering
	filtered_vcf = os.path.join(output_dir, "final_variants/filtered_variants.vcf.gz")
	run_command(f"gatk VariantFiltration -R {reference_genome} -V {raw_vcf} --filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0' --filter-name 'my_snp_filter' -O {filtered_vcf}")

	# Step 10: Variant Annotation with ANNOVAR
	annovar_dir = '/mnt/mone/PMI/CH/02.Variant_Calling/shared/dir-annova/annovar'
	filtered_avinput = os.path.join(output_dir, "final_variants/filtered_variants.avinput")
	annotated_vcf = os.path.join(output_dir, "annotated/filtered_variants_annotated.vcf")

	run_command(f"{annovar_dir}/convert2annovar.pl -format vcf4 {filtered_vcf} -outfile {filtered_avinput}")
	run_command(f"{annovar_dir}/table_annovar.pl {filtered_avinput} {annovar_dir}/CH_hg19_db/ -buildver hg19 -out {output_dir}/annotated/filtered_variants -remove -protocol refGene,cytoBand,gnomad,exac03,avsnp150 -operation gx,r,f,f,f -nastring . -vcfinput")

	print("Variant discovery pipeline completed successfully.")

if __name__=="__main__":
	if len(sys.argv) != 7:
		print("Usage: python variant_call_pipeline.py <input_directory> <output_directory> <reference_genome> <known_sites_vcfs> <total_threads> <samples_per_run>")
		sys.exit(1)

	inputdir = sys.argv[1]
	outputdir = sys.argv[2]
	reference_genome = sys.argv[3]
	known_sites_vcfs = sys.argv[4]
	total_threads = int(sys.argv[5])
	samples_per_run = int(sys.argv[6])
	main(inputdir, outputdir, reference_genome, known_sites_vcfs, total_threads, samples_per_run)


