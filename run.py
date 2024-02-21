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

def domain_files(directory, output_file, output_dir):
    list_file_dict = {}
    bed_file_dict = {}
    conflict_list_dict = {}
    conflict_bed_dict = {}
    species_count_dict = {}
    for filename in os.listdir(directory):
        if filename.endswith(".hmm_results"):
            domain_file_path = os.path.join(directory, filename)
            with open(domain_file_path, 'r') as domain_file:
                for lines in domain_file:
                    if lines.startswith('#') or lines.startswith('-'):
                        continue
                    else:
                        tabbed_lines = re.sub(r'\s+', '\t', lines)
                        columns = tabbed_lines.strip().split('\t')
                        id, domain_name, ali_from, ali_to = columns[0], columns[3], columns[17], columns[18]
                        # Create Key, Value pairs for 'list_file_dict'
                        if domain_name not in list_file_dict:
                            list_file_dict[domain_name] = []
                        list_file_dict[domain_name].append(f'{id}')
                        
                        # Create Key, Value pairs for 'bed_file_dict'
                        if domain_name not in bed_file_dict:
                            bed_file_dict[domain_name] = []
                        bed_file_dict[domain_name].append(f'{id}\t{int(ali_from) - 1}\t{ali_to}')

                        # Create Key, Value pairs for '_Conflict.list'
                        if id not in conflict_list_dict:
                            conflict_list_dict[id] = []
                        conflict_list_dict[id].append(f'{domain_name}')

                        # Create Key Value pairs for '_Conflict.bed'
                        if id not in conflict_bed_dict:
                            conflict_bed_dict[id] = []
                        #conflict_bed_dict[id].append(f'{domain_name}\t{int(ali_from) - 1}\t{ali_to}')
                        conflict_bed_dict[id].append(f'{id}\t{int(ali_from) - 1}\t{ali_to}')
                        
                        # Create a Key Value pairs for table file
                        speciesname = os.path.splitext(filename)[0]
                        if speciesname not in species_count_dict:
                            species_count_dict[speciesname] = {}
                        if domain_name not in species_count_dict[speciesname]:
                            species_count_dict[speciesname][domain_name] = 1
                        else:
                            species_count_dict[speciesname][domain_name] += 1

    for species_name in os.listdir(directory):
        if species_name.endswith(".hmm_results"):
            file_name = species_name.strip().split(".hmm_results")[0]
            
            # Write to list file
            for k, v in list_file_dict.items():
                if output_dir:
                    if not os.path.exists(output_dir):
                        os.makedirs(output_dir)
                    if not output_file:
                        output_list_file = f'{file_name}_{k}.list'
                    else:
                        output_list_file = output_file
                    output_list_file_path = os.path.join(output_dir, output_list_file)
                else:
                    output_list_file_path = output_file or f'{file_name}_{k}.list'
                with open(output_list_file_path, 'w') as list_file:
                    for values in v:
                        list_file.write(f'{values}\n')
    
            # Write to bed file
            for k, v in bed_file_dict.items():
                if output_dir:
                    if not os.path.exists(output_dir):
                        os.makedirs(output_dir)
                    if not output_file:
                        output_bed_file = f'{file_name}_{k}_domain.bed'
                    else:
                        output_bed_file = output_file
                    output_bed_file_path = os.path.join(output_dir, output_bed_file)
                else:
                    output_bed_file_path = output_file or f'{file_name}_{k}_domain.bed'
                with open(output_bed_file_path, 'w') as bed_file:
                    for values in v:
                        bed_file.write(f'{values}\n')

            # Write to Conflict list file
            if output_dir:
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                if not output_file:
                    output_conflict_list_file = f'{file_name}_conflict.list'
                else:
                    output_conflict_list_file = output_file
                output_conflict_list_path = os.path.join(output_dir, output_conflict_list_file)
            else:
                output_conflict_list_path = output_file or f'{file_name}_conflict.list'
            with open(output_conflict_list_path, 'w') as conflict_list:
                for k, v in conflict_list_dict.items():
                    remove_duplicate = list(set(v))
                    if len(remove_duplicate) > 1:
                        #conflicting_domain_list = '\t'.join(remove_duplicate)
                        conflict_list.write(f'{k}\n')

            # Write to Conflict bed file
            if output_dir:
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                if not output_file:
                    output_conflict_bed_file = f'{file_name}_conflict_domain.bed'
                else:
                    output_conflict_bed_file = output_file
                output_conflict_bed_path = os.path.join(output_dir, output_conflict_bed_file)
            else:
                output_conflict_bed_path = output_file or f'{file_name}_conflict_domain.bed'

            with open(output_conflict_bed_path, 'w') as conflict_bed:
                for k, v in conflict_bed_dict.items():
                    if len(v) > 1:
                        for item in v:
                            conflict_bed.write('\t'.join(item.split('\t')) + '\n')
                    #remove_duplicate = list(set(v))
                    #if len(remove_duplicate) > 1:
                        #conflicting_domain_list = '\t'.join(remove_duplicate)
                        #conflict_bed.write(f'{k}\t{conflicting_domain_list}\n')

            # Write to table file
            if output_dir:
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                if not output_file:
                    output_table_file = f'{file_name}_table.txt'
                else:
                    output_table_file = output_file
                output_table_path = os.path.join(output_dir, output_table_file)
            else:
                output_table_path = output_file or f'{file_name}_table.txt'

            sorting_domain_counts = set()
            for name_domain in species_count_dict.values():
                sorting_domain_counts.update(name_domain)

            sorted_domain = sorted(sorting_domain_counts)
            with open(output_table_path, 'w') as table_file:
                table_file.write(f'Species\t')
                for sorted_domain_names in sorted_domain:
                    table_file.write(f'{sorted_domain_names}\t')
                table_file.write(f'\n')
                for k, v in species_count_dict.items():
                    table_file.write(f'{k}\t')
                    for counts in sorting_domain_counts:
                        total = v.get(counts, 0)
                        table_file.write(f'{total}\t')
                    table_file.write(f'\n')

def fasta(directory, secondary, output_file, output_dir):
    for filename in os.listdir(directory):
        if filename.endswith('.list') and not filename.endswith('_domain.list'):
            list_file_path = os.path.join(directory, filename)
            species_name = filename.replace('.list', '')

            for clean_files in os.listdir(secondary):
                if clean_files.endswith('_clean.fasta') and clean_files.startswith(species_name.split('_')[0]):
                    clean_fasta_path = os.path.join(secondary, clean_files)
                    if output_dir:
                        if not os.path.exists(output_dir):
                            os.makedirs(output_dir)
                        if not output_file:
                            output_file_name = f'{species_name}.fasta'
                        else:
                            output_file_name = output_file
                        output_file_name_path = os.path.join(output_dir, output_file_name)
                    else:
                        output_file_name_path = output_file or f'{species_name}.fasta'

                    running_seqkit = f'seqkit grep -f {list_file_path} {clean_fasta_path} -o {output_file_name_path}'
                    subprocess.run(running_seqkit, shell = True, check = True)
        
        elif filename.endswith('_domain.bed'):
            list_file_path = os.path.join(directory, filename)
            species_name = filename.replace('_domain.bed', '')

            for clean_files in os.listdir(secondary):
                if clean_files.endswith('_clean.fasta') and clean_files.startswith(species_name.split('_')[0]):
                    clean_fasta_path = os.path.join(secondary, clean_files)
                    if output_dir:
                        if not os.path.exists(output_dir):
                            os.makedirs(output_dir)
                        if not output_file:
                            output_file_name = f'{species_name}_domain.fasta'
                        else:
                            output_file_name = output_file
                        output_file_name_path = os.path.join(output_dir, output_file_name)
                    else:
                        output_file_name_path = output_file or f'{species_name}_domain.fasta'

                    running_seqkit_subseq = f'seqkit subseq --update-faidx --bed {list_file_path} {clean_fasta_path} -o {output_file_name_path}'
                    subprocess.run(running_seqkit_subseq, shell = True, check = True)

def muscle(directory, output_file, output_dir):
    for domain_fasta in os.listdir(directory):
        if domain_fasta.endswith('_domain.fasta'):
            domain_fasta_path = os.path.join(directory, domain_fasta)
            muscled_name = domain_fasta.replace('_domain.fasta', '_domain_muscled.fasta')
            if output_dir:
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                if not output_file:
                    output_muscle_file_name = f'{muscled_name}'
                else:
                    output_muscle_file_name = output_file
                output_muscle_file_name_path = os.path.join(output_dir, output_muscle_file_name)
            else:
                output_muscle_file_name_path = output_file or f'{muscled_name}'
            
            running_muscle = f'muscle -in {domain_fasta_path} -out {output_muscle_file_name_path}'
            subprocess.run(running_muscle, shell = True, check = True)

def trimal(directory, output_file, output_dir):
    for muscled_fasta_files in os.listdir(directory):
        if muscled_fasta_files.endswith('_domain_muscled.fasta'):
            muscled_files_path = os.path.join(directory, muscled_fasta_files)
            if output_dir:
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                if not output_file:
                    output_trimal_file_name = f'trimal_{muscled_fasta_files}'
                else:
                    output_trimal_file_name = output_file
                output_trimal_file_name_path = os.path.join(output_dir, output_trimal_file_name)
            else:
                output_trimal_file_name_path = output_file or f'trimal_{muscled_fasta_files}'

            running_trimal = f'trimal -in {muscled_files_path} -out {output_trimal_file_name_path} -gappyout'
            subprocess.run(running_trimal, shell = True, check = True)

def weblogo(directory, output_file, output_dir):
    for trimal_and_muscle_files in os.listdir(directory):
        if trimal_and_muscle_files.endswith('_muscled.fasta'):
            trimal_and_muscle_path = os.path.join(directory, trimal_and_muscle_files)
            weblogo_name = trimal_and_muscle_files.replace('_muscled.fasta', '_weblogo.pdf')
            if output_dir:
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                if not output_file:
                    output_trimal_file_name = f'{weblogo_name}'
                else:
                    output_trimal_file_name = output_file
                output_trimal_file_name_path = os.path.join(output_dir, output_trimal_file_name)
            else:
                output_trimal_file_name_path = output_file or f'trimal_{weblogo_name}'
            
            running_weblogo = f'weblogo -f {trimal_and_muscle_path} -D fasta -o {output_trimal_file_name_path} -F pdf --resolution 400'
            subprocess.run(running_weblogo, shell = True, check = True)

def fix_format(directory, output_file, output_dir):
    for fasta_files in os.listdir(directory):
        if fasta_files.endswith('_domain.fasta'):
            fasta_files_path = os.path.join(directory, fasta_files)
            remove_species_name = re.sub(r'^[^_]*_', '', fasta_files_path)
            domain_name = remove_species_name.replace('_domain.fasta', '')
            if output_dir:
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                if not output_file:
                    output_fixed_file_name = f'fixed_{fasta_files}'
                else:
                    output_fixed_file_name = output_file
                output_fixed_file_name_path = os.path.join(output_dir, output_fixed_file_name)
            else:
                output_fixed_file_name_path = output_file or f'fixed_{fasta_files}'

            with open(fasta_files_path, 'r') as reading_fasta, open(output_fixed_file_name_path, 'w') as writing_fixed_fasta:
                for sequence_id in reading_fasta:
                    if sequence_id.startswith('>'):
                        new_header = re.sub(r':.', f'-{domain_name}', sequence_id)
                        writing_fixed_fasta.write(new_header)
                    else:
                        writing_fixed_fasta.write(sequence_id)

def combine_domains(directory, output_file, output_dir):
    species_dict = {}
    for species_file in os.listdir(directory):
        if species_file.startswith('fixed_') and species_file.endswith('_domain.fasta'):
            species_file_path = os.path.join(directory, species_file)
            species_name = species_file.strip().split('_')[1]
            if species_name not in species_dict:
                species_dict[species_name] = []

            with open(species_file_path, 'r') as fasta_information:
                id = ''
                sequence = ''
                for line in fasta_information:
                    if line.startswith('>'):
                        if sequence:
                            species_dict[species_name].append(f'{id}{sequence}')
                        id = line
                        sequence = ''
                    else:
                        sequence += line
                species_dict[species_name].append(f'{id}{sequence}')

    for k, v in species_dict.items():
        if output_dir:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            if not output_file:
                output_fixed_file_name = f'{k}_combined.fasta'
            else:
                output_fixed_file_name = output_file
            output_fixed_file_name_path = os.path.join(output_dir, output_fixed_file_name)
        else:
            output_fixed_file_name_path = output_file or f'{k}_combined.fasta'

        with open(output_fixed_file_name_path, 'w') as output:
            for seq in v:
                output.write(f'{seq}\n')

def main():
    # Initialize the argument parser
    parser = argparse.ArgumentParser(description="List files from specified directories and set number of CPUs.")
    # Add arguments for input folder, database folder, and number of CPUs
    parser.add_argument("-i", "--input", help="Path to the input folder", required=True)
    parser.add_argument("--files", action = "store_true", help = "Usage: python3 run.py -i [Directory Name] --files. Obtain LIST, DOMAIN_BED, CONFLICT_LIST, DOMAIN_CONFLICT_BED and a Table file", required = False)
    parser.add_argument("--fasta", action = "store_true", help = "Usage: python3 run.py -i [Directory Name] --dir [Secondary Folder] --fasta. Obtain FASTA, DOMAIN_FASTA, CONFLICT_FASTA, and DOMAIN_CONFLICT_FASTA files", required = False)
    parser.add_argument("--muscle", action = "store_true", help = "Usage: python3 run.py -i [Directory Name] --muscle.", required = False)
    parser.add_argument("--trimal", action = "store_true", help = "Usage: python3 run.py -i [Directory Name] --trimal.", required = False)
    parser.add_argument("--weblogo", action = "store_true", help = "Usage: python3 run.py -i [Directory Name] --weblogo.", required = False)
    parser.add_argument("--combine", action = "store_true", help = "Usage: python3 run.py -i [Directory Name] --combine.", required = False)
    parser.add_argument("--fix", action = "store_true", help = "Usage: python3 run.py -i [Directory Name] --fix.", required = False)
    parser.add_argument("--db", help="Path to the database folder", required=False)
    parser.add_argument("--CPU", help="Number of CPUs", default=2, type=int, required=False)
    parser.add_argument("--logging", help="Enable logging", action='store_true')

    # Flags
    parser.add_argument("--dir", help = "Option: Directory name", required = False)
    parser.add_argument("--out_file", help = "Option: Name your output File", required = False)
    parser.add_argument("--out_dir", help = "Option: Create output Directory", required = False)
    
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

    if args.files:
        domain_files(args.input, args.out_file, args.out_dir)

    if args.fasta:
        fasta(args.input, args.dir, args.out_file, args.out_dir)

    if args.muscle:
        muscle(args.input, args.out_file, args.out_dir)

    if args.trimal:
        trimal(args.input, args.out_file, args.out_dir)

    if args.weblogo:
        weblogo(args.input, args.out_file, args.out_dir)
    
    if args.fix:
        fix_format(args.input, args.out_file, args.out_dir)

    if args.combine:
        combine_domains(args.input, args.out_file, args.out_dir)

if __name__ == "__main__":
    main()
