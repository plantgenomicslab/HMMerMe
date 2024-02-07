
import os
import argparse
import subprocess
import logging



def process_fasta_files(directory):
	"""
	Processes each FASTA file in the input folder using awk and seqkit.
	"""
	for filename in os.listdir(directory):
		# Ensure filename ends with .fasta or .fa but not with _clean.fasta
		if (filename.endswith(".fasta") or filename.endswith(".fa")) and not filename.endswith("_clean.fasta"):
			base_name = os.path.splitext(filename)[0]  # Extract base name without extension
			input_path = os.path.join(directory, filename)
			clean_path = os.path.join(directory, f"{base_name}_clean.txt")
			clean_fasta_path = os.path.join(directory, f"{base_name}_clean.fasta")

			# Skip processing if the cleaned file already exists
			if os.path.exists(clean_path) and os.path.exists(clean_fasta_path):
				logging.info(f"Skipped processing as cleaned file exists: {clean_path} and {clean_fasta_path}")
				continue

			# Run awk command
			logging.info(f"Processing: {filename}")

			awk_command = f"awk '{{print $1}}' {input_path} > {clean_path}"
			subprocess.run(awk_command, shell=True, check=True)

			seqkit_command = f"seqkit seq {clean_path} -o {clean_fasta_path}"
			subprocess.run(seqkit_command, shell=True, check=True)
			print(f"Processed: {filename}")

def list_files(directory):
	"""
	List all files in a given directory.

	:param directory: The directory path to list files from.
	"""
	if not os.path.exists(directory):
		logging.error(f"The directory {directory} does not exist.")
		return []

	database_path = os.path.join(directory, 'database')
	database_output = open(database_path, "w")

	for filename in os.listdir(directory):
		if filename.endswith(".hmm"):
			input_path = os.path.join(directory, filename)
			base_name = filename.rsplit('.hmm', 1)[0]  # Remove the .hmm extension
			new_hmm_file = os.path.join(directory,f"{base_name}_new")  # New file name
	# Run hmmconvert command
			logging.info(f"Processing HMM file: {input_path}")
			subprocess.run(['hmmconvert', input_path], stdout=database_output)

	database_output.close()
	subprocess.run(f"hmmpress -f {database_path}", shell=True, check=True)
	return database_path


def run_hmm(directory, cpu_count, database_path):
	"""
	Processes each cleaned FASTA file in the input folder using HMMER to search against a database.
	"""
	output_dir = "output"
	os.makedirs(output_dir, exist_ok=True)  # Create output directory if it doesn't exist

	for filename in os.listdir(directory):
		if filename.endswith("_clean.fasta"):
			base_name = filename.rsplit('_clean.fasta', 1)[0]
			input_path = os.path.join(directory, filename)
			base_output_dir = os.path.join(output_dir, base_name)
			os.makedirs(base_output_dir, exist_ok=True)  # Create a subdirectory for each file's output

			output_file_path = os.path.join(base_output_dir, f"{base_name}.hmm_results")
			logging.info(f"Running HMMER for {filename}")
			hmm_command = ["hmmsearch", "--noali", "--domE", "1e-10", "--domtblout", output_file_path, "-E", "1e-10", "--cpu", str(cpu_count), database_path, input_path]
#			subprocess.run(hmm_command , check=True)
			
			try:
				result = subprocess.run(hmm_command, check=True, capture_output=True, text=True)
				logging.info(result.stdout)
			except subprocess.CalledProcessError as e:
				logging.error(f"Error running HMMER: {e}")
				logging.error(e.stderr)

####
# for i in `ls *hmm | sed 's/.hmm//g'`;
#   do hmmconvert ${i}.hmm > ${i}_new.hmm;
# done
# cat *_new.hmm >> database.hmm

# hmmpress database.hmm
###

def main():
	# Initialize the argument parser
	parser = argparse.ArgumentParser(description="List files from specified directories and set number of CPUs.")
	# Add arguments for input folder, database folder, and number of CPUs
	parser.add_argument("--input", help="Path to the input folder", required=True)
	parser.add_argument("--db", help="Path to the database folder", required=True)
	parser.add_argument("--CPU", help="Number of CPUs", type=int, required=True)
	parser.add_argument("--logging", help="Enable logging", action='store_true')

	# Parse the arguments
	args = parser.parse_args()

	if args.logging:
		logging.basicConfig(level=logging.INFO,filename='hmm_run.log', format='%(asctime)s - %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
	else:
		logging.basicConfig(level=logging.CRITICAL)  # Print CRI
		
	logging.info("Starting processing of FASTA files.")
	process_fasta_files(args.input)

	logging.info("Creating HMM database.")
	database_path = list_files(args.db)
	if database_path:
		logging.info("Running HMM searches.")
		run_hmm(args.input, args.CPU, database_path)
	else:
		logging.error("Database creation failed or was skipped.")

if __name__ == "__main__":
	main()
