#!/usr/bin/env python
import os, argparse, subprocess, logging, re

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
    if directory is None:
        logging.error("No database directory provided.")
        return None
    
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
            # Run hmmsearch command and log the result
            try:
                result = subprocess.run(hmm_command, check=True, capture_output=True, text=True)
                logging.info(result.stdout)
            except subprocess.CalledProcessError as e:
                logging.error(f"Error running HMMER: {e}")
                logging.error(e.stderr)

def domain_files(input, output=None):

    # Non - conflict
    list_file_dict = {}
    bed_file_dict = {}

    # Conflict
    conflict_list_dict = {}
    conflict_bed_dict = {}

    # Table
    species_count_dict = {}

    # Dictonary function
    for hmm_result_file in os.listdir(input):
        if hmm_result_file.endswith('.hmm_results'):
            hmm_result_path = os.path.join(input, hmm_result_file)

            # Assign outputfile name
            name = hmm_result_file.strip().split('.hmm_results')[0]

            with open(hmm_result_path, 'r') as is_hmm_file:
                for lines in is_hmm_file:
                    if lines.startswith('#') or lines.startswith('-'):
                        continue
                    else:
                        # Assign variables
                        species_name = os.path.splitext(hmm_result_file)[0]
                        remove_blanks_from_lines_and_replace_with_tab = re.sub(r'\s+', '\t', lines)
                        columns = remove_blanks_from_lines_and_replace_with_tab.strip().split('\t')
                        gene_id, domain_id, ali_from, ali_to = columns[0], columns[3], columns[17], columns[18]

                        # Creating list_file_dict: Key = domain_id, Value = gene_id
                        if domain_id not in list_file_dict:
                            list_file_dict[domain_id] = []
                        list_file_dict[domain_id].append(gene_id)

                        # Creating bed_file_dict: Key = domain_id, Value = gene_id'\t'ali_from'\t'ali_to
                        if domain_id not in bed_file_dict:
                            bed_file_dict[domain_id] = []
                        bed_file_dict[domain_id].append(f'{gene_id}\t{ali_from}\t{ali_to}')

                        # Creating conflict_list_dict: Key = gene_id, Value = domain_id
                        if gene_id not in conflict_list_dict:
                            conflict_list_dict[gene_id] = []
                        conflict_list_dict[gene_id].append(domain_id)

                        # Creating conflict_bed_dict: Key = gene_id, Value = domain_id
                        if gene_id not in conflict_bed_dict:
                            conflict_bed_dict[gene_id] = []
                        conflict_bed_dict[gene_id].append(f'{gene_id}\t{ali_from}\t{ali_to}')

                        # Creating species_count_dict: Key = species_name[domain_id], Value = {}
                        if species_name not in species_count_dict:
                            species_count_dict[species_name] = {}
                        if domain_id not in species_count_dict[species_name]:
                            species_count_dict[species_name][domain_id] = 1
                        else:
                            species_count_dict[species_name][domain_id] += 1
    
            # '.list' file, directory
            for k, v in list_file_dict.items():
                list_output_name = f'{name}_{k}.list'
                if output:
                    if not os.path.exists(output):
                        os.makedirs(output)
                else:
                    pass
                with open(os.path.join(output, list_output_name), 'w') as list_file:
                    for gene_id in v:
                        list_file.write(f'{gene_id}\n')

            # '.bed' file, directory
            for k, v in bed_file_dict.items():
                bed_output_name = f'{name}_{k}.bed'
                if output:
                    if not os.path.exists(output):
                        os.makedirs(output)
                else:
                    pass
                with open(os.path.join(output, bed_output_name), 'w') as bed_file:
                    for geneID_aliFrom_aliTo in v:
                        bed_file.write(f'{geneID_aliFrom_aliTo}\n')
            
            # '_conflict.list' file, directory
            conflict_list_name = f'{name}_conflict.list'
            if output:
                if not os.path.exists(output):
                    os.makedirs(output)
            else:
                pass
            with open(os.path.join(output, conflict_list_name), 'w') as conflict_list_file:
                for k, v in conflict_list_dict.items():
                    remove_duplicates_domain_id = list(set(v))
                    if len(remove_duplicates_domain_id) > 1:
                        conflict_list_file.write(f'{k}')

            # '_conflict_domain.bed' file, directory
            conflict_domain_bed_name = f'{name}_conflict_domain.bed'
            if output:
                if not os.path.exists(output):
                    os.makedirs(output)
            else:
                pass
            with open(os.path.join(output, conflict_domain_bed_name), 'w') as conflict_bed_file:
                for k, v in conflict_bed_dict.items():
                    if len(v) > 1:
                        for bed_format in v:
                            conflict_bed_file.write('\t'.join(bed_format.split('\t')) + '\n')

def main():
    # Initialize the argument parser
    parser = argparse.ArgumentParser(description="List files from specified directories and set number of CPUs.")
    # Add arguments for input folder, database folder, and number of CPUs
    parser.add_argument("--input", help="Path to the input folder", required=True)
    parser.add_argument("--output", help="Create Output directory, default is 'output'.", default='.', required=False)
    parser.add_argument("--db", help="Path to the database folder", required=False)
    parser.add_argument("--CPU", help="Number of CPUs", type=int, required=False)
    parser.add_argument("--logging", help="Enable logging", action='store_true')

    # Flags1
    parser.add_argument("--weblogo", help="Create Weblogo Design", action='store_true', required=False)

    # Parse the arguments
    args = parser.parse_args()

    # Configure logging based on user input
    if args.logging:
        logging.basicConfig(level=logging.INFO,filename='hmm_run.log', format='%(asctime)s - %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    else:
        logging.basicConfig(level=logging.CRITICAL)  # Print CRITICAL logs only

    logging.info("Starting processing of FASTA files.")
    process_fasta_files(args.input)

    logging.info("Creating HMM database.")
    database_path = list_files(args.db)
    if database_path:
        logging.info("Running HMM searches.")
        run_hmm(args.input, args.CPU, database_path)
    else:
        logging.error("Database creation failed or was skipped.")

    domain_files(args.input, args.output)

if __name__ == "__main__":
    main()
