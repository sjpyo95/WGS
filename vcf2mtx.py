import sys
import csv

def convert_gt(gt):
	if gt == './.': return ''
	elif gt == '0/0': return '0'
	elif gt in ('1/0', '0/1'): return '1'
	elif gt == '1/1': return '2'
	else: return ''

def vcf2mtx(vcf, outfile):
	rows = []
#	sampids = []

	with open(vcf, 'r') as f:
		for line in f:
			line = line.strip()
			if line.startswith('##'): continue
			if line.startswith('#CHROM'):
				header = line.split('\t')
				sampids = header[9:]
				continue
			tmp = line.split('\t')
			chrom = tmp[0]
			pos = tmp[1]
			ref = tmp[3]
			alt = tmp[4]
			varkey = f'{chrom}:{pos}|{ref}->{alt}'
			sampvals = []
			for sampval in tmp[9:]:
				gtstr = sampval.split(':')[0]
				sampvals.append(convert_gt(gtstr))
			row = [varkey] + sampvals
			rows.append(row)
	with open(outfile, 'w', newline='') as csvfile:
		writer = csv.writer(csvfile)
		writer.writerow(['variant']+sampids)
		writer.writerows(rows)
	print(f'matrix CSV saved to {outfile}')

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print('Usage: python vcf2mtx.py <merged_vcf_file> <output_csv>')
		sys.exit(1)
	vcf = sys.argv[1]
	outfile = sys.argv[2]
	vcf2mtx(vcf, outfile)
