# WGS
Protocol for Finding and Organizing Germline Variants Correlated with CHIP and Measuring CH Risk Using Machine Learning

## Step 1: Data Preprocessing
### 1. Download and Prepare WGS Data:
		○ Access the WGS data from the specified project accession (PRJEB50398).
		○ Ensure all samples are correctly formatted and paired-end sequences are identified.		

###	2. Quality Control (QC):
		○ Use tools like FastQC to check the quality of the sequencing reads.
		○ Trim low-quality bases and adapters using tools like Trimmomatic or Cutadapt.

###	3. Align Reads to Reference Genome:
		○ Use BWA or a similar aligner to align the reads to the human reference genome (e.g., GRCh38).
###	4. Mark Duplicates and Recalibrate Base Quality Scores:
		○ Use tools like Picard to mark duplicates.

## Step 2: Variant Calling and Annotation
### 1. Call Germline Variants:
		○ Use GATK HaplotypeCaller to call germline variants from the aligned reads.
		○ Perform joint genotyping if necessary.
	
 ### 2. Variant Filtering:
		○ Apply variant quality score recalibration (VQSR) using GATK or hard filtering criteria to filter out low-quality variants.
	
 ### 3. Annotate Variants:
Annotate the variants using tools like ANNOVAR, SnpEff, or VEP to obtain functional information about the variants.
