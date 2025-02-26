import os
import subprocess
import zipfile
import argparse
from pathlib import Path

def extract_fastqc_data(fastqc_zip_file, outputdir):
	"""
	Extract fastqc_data.txt from the FastQC zip file.
	"""
	with zipfile.ZipFile(fastqc_zip_file, 'r') as zip_ref:
		# Extract the fastqc_data.txt file from the ZIP
		for file in zip_ref.namelist():
			if "fastqc_data.txt" in file:
				zip_ref.extract(file, outputdir)
				return os.path.join(outputdir, file)
	raise FileNotFoundError("fastqc_data.txt not found in the ZIP file.")

def parse_fastqc(fastqc_data_file):
	"""
	Parse FastQC result file to extract Q20, Q30, total reads, and read length.
	"""
	fastqc_data = {
		"% > Fastq Q20": None,
		"% > Fastq Q30": None,
		"Total reads": None,
		"Read length (bp)": None,
	}
	
	total_seq = 0
	q20_seq = 0
	q30_seq = 0

	with open(fastqc_data_file, 'r') as f:
		in_q_section = False

		for line in f:
			if "Total Sequences" in line:
				fastqc_data["Total reads"] = int(line.strip().split('\t')[1])
			elif "Sequence length" in line:
				fastqc_data["Read length (bp)"] = int(line.strip().split('\t')[1])
			elif ">>Per sequence quality scores" in line:
				in_q_section = True
				continue
			elif ">>END_MODULE" in line and in_q_section:
				in_q_section = False

			if in_q_section and not line.startswith("#"):
				parts = line.strip().split('\t')
				quality = float(parts[0])
				counts = float(parts[1])
				
				total_seq += counts
				if quality >= 20:
					q20_seq += counts
				if quality >= 30:
					q30_seq += counts


		if total_seq > 0:
			fastqc_data['% > Fastq Q20'] = (q20_seq / total_seq) * 100
			fastqc_data['% > Fastq Q30'] = (q30_seq / total_seq) * 100

		return fastqc_data

def run_samtools_stats(bam_file):
	"""
	Run samtools to extract BAM file statistics.
	"""
	samtools_stats_cmd = f"samtools stats {bam_file}"
	result = subprocess.run(samtools_stats_cmd, shell=True, check=True, capture_output=True, text=True)

	stats_data = {
		"Total reads": None,
		"Mappable reads": None,
		"Mappable yield (Mbp)": None,
		"Throughput mean depth (X)": None,
		"Fragment length median": None,
		"Standard deviation of insert length": None
	}

	for line in result.stdout.splitlines():
		if 'raw total sequences:' in line:
			stats_data['Total reads'] = int(line.split(':')[1].strip())
		elif 'reads mapped:' in line:
			stats_data['Mappable reads'] = int(line.split(':')[1].strip())
		elif 'insert size average' in line:
			stats_data['Fragment length median'] = float(line.split(':')[1].strip())
		elif 'insert size standard deviation' in line:
			stats_data['Standard deviation of insert length'] = float(line.split(':')[1].strip())

	return stats_data

def run_picard_dedup(bam_file, outputdir):
	"""
	Run Picard MarkDuplicates to get de-duplication metrics.
	"""
	os.makedirs(outputdir, exist_ok=True)
	dedup_metrics_file = os.path.join(outputdir, "dedup_metrics.txt")
	dedup_bam = os.path.join(outputdir, "dedup.bam")

	picard_cmd = f"picard MarkDuplicates I={bam_file} O={dedup_bam} M={dedup_metrics_file} REMOVE_DUPLICATES=true"
	subprocess.run(picard_cmd, shell=True, check=True)

	dedup_data = {
		"De-duplicated reads": None,
		"De-duplicated reads % (out of total reads)": None,
	}

	with open(dedup_metrics_file, 'r') as f:
		for line in f:
			if "READS_UNIQUE" in line:
				dedup_data["De-duplicated reads"] = int(line.split('\t')[2].strip())
			elif 'PERCENT_DUPLICATION' in line:
				dedup_data['De-duplicated reads % (out of total reads)'] = str((1 - float(line.split("\t")[1].strip())) * 100)

	return dedup_data

def run_bedtools_coverage(bam_file, ref_genome):
	"""
	Run bedtools to calculate coverage at different levels (1X, 5X, 10X, etc.)
	"""
	coverage_data = {}
	genome_size_cmd = f"samtools faidx {ref_genome} | awk '{{sum += $2}} END {{print sum}}'"
	genome_size_result = subprocess.run(genome_size_cmd, shell=True, check=True, capture_output=True, text=True)
	genome_size = int(genome_size_result.stdout.strip())

	for coverage in [1, 5, 10, 15, 20, 30]:
		bedtools_cmd = f"bedtools genomecov -ibam {bam_file} -bga | awk '{{if ($4 >= {coverage}) sum += ($3 - $2)}} END {{print sum}}'"
		result = subprocess.run(bedtools_cmd, shell=True, check=True, capture_output=True, text=True)

		covered_bases = int(result.stdout.strip())
		coverage_percentage = (covered_bases / genome_size) * 100
		coverage_data[f"% >= {coverage}X coverage"] = round(coverage_percentage, 2)

	return coverage_data

def main():
	parser = argparse.ArgumentParser(description="Calculate QC metrics for sequencing data.")
	parser.add_argument("--fastqc_zip", nargs='+', required=True, help="Path to the FastQC ZIP files.")
	parser.add_argument("--bam_file", required=True, help="Path to the sorted BAM file.")
	parser.add_argument("--ref_genome", required=True, help="Path to the reference genome file.")
	parser.add_argument("--outputdir", required=True, help="Directory to save the outputs.")
	args = parser.parse_args()

	outputdir = args.outputdir
	os.makedirs(outputdir, exist_ok=True)

	# Process FastQC data
	fastqc_data_files = [extract_fastqc_data(args.fastqc_zip[0], outputdir), extract_fastqc_data(args.fastqc_zip[1], outputdir)]
	fastqc_data1 = parse_fastqc(fastqc_data_files[0])
	fastqc_data2 = parse_fastqc(fastqc_data_files[1])

	fastqc_data = {
		"% > Fastq Q20": round((fastqc_data1["% > Fastq Q20"] + fastqc_data2["% > Fastq Q20"]) / 2, 1),
		"% > Fastq Q30": round((fastqc_data1["% > Fastq Q30"] + fastqc_data2["% > Fastq Q30"]) / 2, 1),
		"Total reads": fastqc_data1["Total reads"] + fastqc_data2["Total reads"],
		"Read length (bp)": round((fastqc_data1["Read length (bp)"] + fastqc_data2["Read length (bp)"]) / 2, 1),
		}

		
	# Process BAM statistics
	bam_stats = run_samtools_stats(args.bam_file)

	# Process deduplication stats
	dedup_stats = run_picard_dedup(args.bam_file, outputdir)

	# Process coverage stats
	coverage_stats = run_bedtools_coverage(args.bam_file, args.ref_genome)

	# Combine all metrics
	all_metrics = {**fastqc_data, **bam_stats, **dedup_stats, **coverage_stats}

	# Save metrics to a file
	output_file = os.path.join(outputdir, "qc_metrics.csv")
	with open(output_file, 'w', newline='') as csvfile:
		fieldnames = all_metrics.keys()
		writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
		writer.writeheader()
		writer.writerow(all_metrics)

	print(f"QC metrics saved to {output_file}")

if __name__ == "__main__":
	main()

