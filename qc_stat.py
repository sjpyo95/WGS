import sys, os
import pysam
import zipfile
import subprocess
from collections import defaultdict

def parse_fastqc(fastqcZip):
	with zipfile.ZipFile(fastqcZip, 'r') as zipref:
		for file in zipref.namelist():
			if 'fastqc_data.txt' in file:
				datafile = file
		with zipref.open(datafile) as file:
			totalcount = None
			q20n = 0
			q30n = 0
			qualsec = False

			for line in file:
				line = line.decode('utf-8').strip()

				if line.startswith('Total Sequences'):
					totalcount = int(line.split('\t')[1])

				if line.startswith('Sequence length'):
					seqlen = int(line.split('\t')[1].split('-')[1])

				if line.startswith('>>Per sequence quality scores'):
					qualsec = True
					continue

				if line.startswith('>>END_MODULE') and qualsec:
					break

				if qualsec and not line.startswith('#'):
					parts = line.split('\t')
					score = float(parts[0])
					count = float(parts[1])
					
					if score >= 20:
						q20n += count

					if score >= 30:
						q30n += count

			q20per = (q20n / totalcount ) * 100 
			q30per = (q30n / totalcount ) * 100
			
			totalyield = totalcount * seqlen / 1000000	# Chanage to Mbp
			print('reads stat done')
			return q20per, q30per, totalcount, seqlen, totalyield

def cal_refsize(faidxfile):
	total = 0
	with open(faidxfile, 'r') as f:
		for line in f:
			l = int(line.split('\t')[1])
			total += l
	refsize = total / 1000000 # Change to Mbp
	print('refsize done')
	return refsize

def mapping_stat(bamfile):
	flagstatOut = subprocess.check_output(f'samtools flagstat {bamfile}', shell=True).decode('utf-8')
	for line in flagstatOut.splitlines():
		if 'primary' in line:
			totalreads = int(line.split()[0])
		elif 'mapped' in line:
			mappedreads = int(line.split()[0])
		elif 'duplicates' in line:
			dupreads = int(line.split()[0])

	dedupreads = totalreads - dupreads
	dedupreads_per = dedupreads / totalreads * 100

	mappedreads_per = (mappedreads/dedupreads) * 100 if dedupreads else 0
	print('mapping stat done')
	return dedupreads, dedupreads_per, mappedreads, mappedreads_per

def cal_coverage(bamfile, ths=[1,5,10,15,20,30]):
	bam = pysam.AlignmentFile(bamfile, 'rb')

	genome_size = sum(bam.lengths)

	coverage_all = defaultdict(int)
	coverage_dedup = defaultdict(int)
	coverage_unique = defaultdict(int)

	for pileup_col in bam.pileup():
		pos = pileup_col.pos

		depth_all = 0
		depth_dedup = 0
		depth_unique = 0

		for pileup_read in pileup_col.pileups:
			if not pileup_read.is_del and not pileup_read.is_refskip:
				alignment = pileup_read.alignment

				depth_all += 1

				if not alignment.is_duplicate:
					depth_dedup += 1

				if not alignment.is_duplicate and alignment.mapping_quality > 0:
					depth_unique += 1

		coverage_all[depth_all] += 1
		coverage_dedup[depth_dedup] += 1
		coverage_unique[depth_unique] += 1

	results = {}
	for category, coverage_dict in [
		('all', coverage_all),
		('dedup', coverage_dedup),
		('unique', coverage_unique)
	]:
		for th in ths:
			covered_bases = sum(count for depth, count in coverage_dict.items() if depth >= th)
			percentage = (covered_bases / genome_size) * 100
			results[f'{category}_{th}x'] = round(percentage, 2)

	bam.close()
	return results

def cal_fragmentStat(bamfile):
	import statistics
	insertsizes = []

	with pysam.AlignmentFile(bamfile, 'rb') as bam:
		for read in bam:
			if read.is_proper_pair and not read.is_secondary and read.template_length > 0:
				insertsizes.append(read.template_length)
	
	fragMedian = statistics.median(insertsizes)
	fragStdev = statistics.stdev(insertsizes)
	print('fragment stat done')
	return fragMedian, fragStdev

def main():
	fastqc1 = sys.argv[1]
	fastqc2 = sys.argv[2]
	bamfile = sys.argv[3]
	outfile = sys.argv[4]
	
	fastqcZips = [fastqc1, fastqc2]

	tq20p = 0; tq30p = 0; tcount = 0; tseqlen = 0; tyield = 0
	for fastqcZip in fastqcZips:
		q20p, q30p, count, seqlen, subyield = parse_fastqc(fastqcZip)
		tq20p += q20p
		tq30p += q30p
		tcount += count
		tseqlen += seqlen
		tyield += subyield
	
	q20per = tq20p / len(fastqcZips)
	q30per = tq30p / len(fastqcZips)
	readlen = tseqlen / len(fastqcZips)
	
	refsize = cal_refsize('/home/sunyme95/CTCELLS/wgs_pipeline/workflow/reference/bwa-mem/genome.fa.fai')
	tmd = tyield / refsize

	dedupreads, dedupreads_per, mappedreads, mappedreads_per = mapping_stat(bamfile)

	mapyield = (mappedreads * readlen) / 1000000
	
	mmd = mapyield / refsize

	ths=[1,5,10,15,20,30]
	covperDic = cal_coverage(bamfile)
	print(covperDic)

	fragMedian, fragStdev = cal_fragmentStat(bamfile)

	sampleid = os.path.basename(bamfile).split('.dedup')[0]
	with open(outfile, 'w') as f:
		f.write('\n')
		f.write('\t'.join(map(str,[sampleid, round(q20per,2), round(q30per,2), tcount, readlen, round(tyield,2), round(refsize,0), round(tmd, 2),dedupreads, round(dedupreads_per, 2), mappedreads, round(mappedreads_per, 2), round(mapyield, 2), round(mmd), *[round(covperDic[t],2) for t in ths], fragMedian, round(fragStdev, 2)])))

if __name__=="__main__":
	main()
