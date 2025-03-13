import sys, os
import pandas as pd

def parseInfoCsv(infoCsvFile):
	df = pd.read_csv(infoCsvFile)
	infoDict = {}
	for idx, row in df.iterrows():
		sampleID = row['Sample ID']
		vcfName = row['VCF name']
		if '.CEL' in vcfName:
			origName = vcfName.split('.CEL')[0]

		else:
			origName = vcfName
		group = row['Group']
		age = 2025 - int(str(row['Birthday'])[:4])
		gender = row['Gender']
		newSampleName = f"{group}-{age}-{gender}-{sampleID}"
		infoDict[vcfName] = {
			'SampleID': sampleID,
			'OriginalSampleName': origName,
			'NewSampleName': newSampleName,
			'Group': group,
			'Age': age,
			'Gender': gender
		}
	return infoDict, df

def writeSampleInfoTable(infoDict, df, outputPrefix):

	rows = []
	for vcfName, info in infoDict.items():
		rows.append({
			'Sample ID': info['SampleID'],
			'Original Sample Name': info['OriginalSampleName'],
			'New Sample Name': info['NewSampleName'],
			'Group': info['Group'],
			'Age': info['Age'],
			'Gender': info['Gender']
		})
	outDf = pd.DataFrame(rows)
	outFile = f'{outputPrefix}_sampleInfo.csv'
	outDf.to_csv(outFile, index=False)
	print(f'Sample Information table saved to {outFile}')

def recodeGenotype(gt, alleleIndex):
	if not ('/' in gt):
		if gt == '0': return '0/0'
		elif gt == '1': return '1/1'
		else: return './.'
	alleles = gt.split('/')
	new_alleles = []
	for a in alleles:
		if a == '0':
			new_alleles.append('0')
		elif a == str(alleleIndex):
			new_alleles.append('1')
		else:
			new_alleles.append('0')
	return '/'.join(new_alleles)
	
def updateInfoField(info, alleleIndex):
	parts = info.split(';')
	new_parts = []
	for part in parts:
		if part.startswith('AC='):
			val = part[3:]
			ac_vals = val.split(',')
			if len(ac_vals) >= alleleIndex:
				new_ac = ac_vals[alleleIndex - 1]
			else:
				new_ac = ''
			new_parts.append(f'AC={new_ac}')
		else:
			new_parts.append(part)
	return ';'.join(new_parts)

def processVcfLine(line):
	fields = line.strip().split('\t')
	if len(fields) < 10:
		return []
	chrom = fields[0]
	pos = fields[1]
	origID = fields[2]
	ref = fields[3]
	alt_field = fields[4]
	qual = fields[5]
	filt = fields[6]
	info = fields[7]
	fmt =fields[8]
	gt = fields[9].split(':')[0]

	newID = origID
	for item in info.split(';'):
		if item.startswith('RSID='):
			newID = item.split('=')[1]
			break
	alt_alleles = alt_field.split(',')
	entries = []
	for i, allele in enumerate(alt_alleles, start=1):
		new_gt = recodeGenotype(gt, i)
		new_info = updateInfoField(info, i)
		entry = {
			'CHROM': chrom,
			'POS': pos,
			'ID': newID,
			'REF': ref,
			'ALT': allele,
			'QUAL': qual,
			'FILTER': filt,
			'INFO': new_info,
			'FORMAT': fmt,
			'GT': new_gt
		}
		entries.append(entry)
	return entries

def mergeVcf(vcfs, infoDict, outputPrefix):
	merged_variants = {}
	sample_names = []
	for vcfPath in vcfs:
		base = os.path.basename(vcfPath)
		if base not in infoDict:
			print(f'Warning: {base} not found in info dictionary. Skipping.')
			continue
		sample_info = infoDict[base]
		newSampleName = sample_info['NewSampleName']
		sample_names.append(newSampleName)
		with open(vcfPath) as f:
			for line in f:
				if line.startswith('#'):
					continue
				entries = processVcfLine(line)
				for entry in entries:
					varKey = f"{entry['CHROM']}:{entry['POS']}|{entry['REF']}|{entry['ALT']}"
					if varKey not in merged_variants:
						merged_variants[varKey] = {
							'CHROM': entry['CHROM'],
							'POS': entry['POS'],
							'ID': entry['ID'],
							'REF': entry['REF'],
							'ALT': entry['ALT'],
							'QUAL': entry['QUAL'],
							'FILTER': entry['FILTER'],
							'INFO': entry['INFO'],
							'FORMAT': entry['FORMAT'],
							'genotypes': {}
						}
					merged_variants[varKey]['genotypes'][newSampleName] = entry['GT']

	def sort_key(key):
		chrom, rest = key.split(':', 1)
		pos = int(rest.split('|')[0])
		return (chrom, pos)
	sorted_keys = sorted(merged_variants.keys(), key=sort_key)

	merged_vcf_file = f"{outputPrefix}.vcf"
	with open(merged_vcf_file, 'w') as out:
		out.write("##fileformat=VCFv4.2\n")
		out.write("##source=MergedVCF\n")
		header_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
		for s in sample_names:
			header_line += f"\t{s}"
		header_line += '\n'
		out.write(header_line)

		for key in sorted_keys:
			var = merged_variants[key]
			ac = 0
			ns = 0

			for s in sample_names:
				gt = var['genotypes'].get(s, './.')
				if gt == './.':
					continue
				ns += 1

				if gt in ('0/1', '1/0'):
					ac += 1
				elif gt == '1/1':
					ac += 2
			an = ns * 2
			var['INFO'] = f"AC={ac};AN={an};NS={ns}"

			genotype_list = []
			for s in sample_names:
				gt = var['genotypes'].get(s, './.')
				genotype_list.append(gt)
			line_fields = [
				var['CHROM'],
				var['POS'],
				var['ID'],
				var['REF'],
				var['ALT'],
				var['QUAL'],
				var['FILTER'],
				var['INFO'],
				var['FORMAT']
			] + genotype_list
			out.write('\t'.join(line_fields) + '\n')
	print(f'Merged VCF file saved to {merged_vcf_file}')

def main():
	if len(sys.argv) != 4:
		print('Usage: python mergeVcf.py <vcfDir> <infoCsvFile> <outputPrefix>')
		sys.exit(1)
	vcfDir = sys.argv[1]
	infoCsvFile = sys.argv[2]
	outputPrefix = sys.argv[3]

	infoDict, infoDf = parseInfoCsv(infoCsvFile)
	writeSampleInfoTable(infoDict, infoDf, outputPrefix)

	vcfs = []
	for f in os.listdir(vcfDir):
		if f.endswith('.vcf'):
			vcfs.append(os.path.join(vcfDir, f))
	if not vcfs:
		print('No VCF files found in the directory.')
		sys.exit(1)
	
	mergeVcf(vcfs, infoDict, outputPrefix)

if __name__ == '__main__':
	main()
