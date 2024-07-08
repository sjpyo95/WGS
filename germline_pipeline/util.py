import subprocess
import pandas as pd
import os
from glob import glob
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

def run_command(command):
	"""Utility function to run a command in the shell."""
	process = subprocess.run(command, shell=True, check=True, text=True, capture_output=True)
	print(process.stdout)
	if process.stderr:
		print(process.stderr)

def read_table(infile):
	"""Read the sample information file and return a DataFrame."""
	return pd.read_csv(infile, sep='\t')

def create_multi_dir(output_dir, directories):
	"""Create necessary output directories."""
	for directory in directories:
		subprocess.run(f"mkdir -p {output_dir}/{directory}", shell=True)

def run_parallel(sample_dirs, input_dir, output_dir, total_threads, samples_per_run, process_func, *args):
	"""Run the pipeline in parallel for multiple samples."""
	threads_per_sample = total_threads // samples_per_run

	total_samples = len(sample_dirs)
	completed_samples = 0

	results = []
	with ThreadPoolExecutor(max_workers=samples_per_run) as executor:
		futures = {executor.submit(process_func, sample_id, input_dir, output_dir, threads_per_sample, *args): sample_id for sample_id in sample_dirs}
		
		for future in tqdm(as_completed(futures), total=len(futures), unit="sample", ncols=100, desc="Processing samples", bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} {percentage:3.0f}%"):
			sample_id = futures[future]
			try:
				result = future.result()
				if result:
					completed_samples += 1
					results.append(result)
				print(f"Sample {sample_id} processed successfully.")
			except Exception as e:
				print(f"Sample {sample_id} failed with error: {e}")

			progress_percentage = (completed_samples / total_samples) * 100
			print(f"Completed {completed_samples}/{total_samples} ({progress_percentage:.2f}%)")

	return results
