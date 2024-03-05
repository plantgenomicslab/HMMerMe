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

        # Using subprocess.Popen for commands that involve redirections
            with open(clean_path, "w") as outfile:
                subprocess.run(["awk", "{print $1}", input_path], stdout=outfile, check=True)


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

def run_hmm(directory, cpu_count, database_path, weblogo):
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

    #Output/Species_folder/hmm_results
    for species_folder_in_output in os.listdir(output_dir):
        species_folder_path = os.path.join(output_dir, species_folder_in_output)

        if not os.path.isdir(species_folder_path):
            continue

        # Non conflict
        list_file_dict = {}
        bed_file_dict = {}    

        # Conflict
        conflict_list_dict = {}
        conflict_bed_dict = {}

        # Table
        table_count_dict = {}

        for hmm_result_file in os.listdir(species_folder_path):
            if hmm_result_file.endswith('.hmm_results'):
                species_name = hmm_result_file.strip().split('.hmm_results')[0]
                hmm_result_file_path = os.path.join(species_folder_path, hmm_result_file)
                with open(hmm_result_file_path, 'r') as result_file:
                    for lines in result_file:
                        if lines.startswith('#') or lines.startswith('-'):
                            continue
                        else:
                            remove_empty_spaces = re.sub(r'\s+', '\t', lines)
                            columns = remove_empty_spaces.strip().split('\t')
                            gene_id, domain_id, ali_from, ali_to = columns[0], columns[3], columns[17], columns[18]
                            
                            # Create list_file_dict, Keys = domain_id : Values = gene_id
                            if domain_id not in list_file_dict:
                                list_file_dict[domain_id] = []
                            list_file_dict[domain_id].append(gene_id)

                            # Create bed_file_dict, Keys = domain_id : Values = gene_id \t ali_from \t ali_to
                            if domain_id not in bed_file_dict:
                                bed_file_dict[domain_id] = []
                            bed_file_dict[domain_id].append(f'{gene_id}\t{int(ali_from) - 1}\t{ali_to}')

                            # Create conflict_list_dict, keys = gene_id : values = domain_id
                            if gene_id not in conflict_list_dict:
                                conflict_list_dict[gene_id] = []
                            conflict_list_dict[gene_id].append(domain_id)

                            # Create conflict_bed_dict, Keys = gene_id : values = ali_from \t ali_to
                            if gene_id not in conflict_bed_dict:
                                conflict_bed_dict[gene_id] = []
                            conflict_bed_dict[gene_id].append(f'{int(ali_from) - 1}\t{ali_to}') 

                            # Create Table, Keys = species_name : Values = counts
                            if species_name not in table_count_dict:
                                table_count_dict[species_name] = {}
                            if domain_id not in table_count_dict[species_name]:
                                table_count_dict[species_name][domain_id] = 1
                            else:
                                table_count_dict[species_name][domain_id] += 1

                # Write '{species}_{domain}.list' files
                for k, v in list_file_dict.items():
                    output_list_file_name = f'{species_name}_{k}.list'
                    print('-' * len(output_list_file_name) * 5)
                    print(f'Preparing to Write in: {output_list_file_name}')
                    print('\n')
                    with open(os.path.join(species_folder_path, output_list_file_name), 'w') as writing_list_output_file:
                        for genes in v:
                            writing_list_output_file.write(f'{genes}\n')
                        print(f'Successfully written your {output_list_file_name} in directory: {species_folder_in_output}')
                        print('-' * len(output_list_file_name) * 5)
                        print('\n')

                # Write '{species}_{domain}_domain.bed' files
                for k, v in bed_file_dict.items():
                    output_bed_file_name = f'{species_name}_{k}_domain.bed'
                    print('-' * len(output_bed_file_name) * 5)
                    print(f'Preparing to write in: {output_bed_file_name}')
                    print('\n')
                    with open(os.path.join(species_folder_path, output_bed_file_name), 'w') as writing_bed_output_file:
                        for bed_format in v:
                            writing_bed_output_file.write(f'{bed_format}\n')
                        print(f'Successfully written your {output_bed_file_name} in directory: {species_folder_in_output}')
                        print('-' * len(output_bed_file_name) * 5)
                        print('\n')

                # Write '{species}_conflict.list' files
                output_conflict_file_name = f'{species_name}_conflict.list'
                print('-' * len(output_conflict_file_name) * 5)
                print(f'Preparing to write in: {output_conflict_file_name}')
                print('\n')
                with open(os.path.join(species_folder_path, output_conflict_file_name), 'w') as writing_conflict_output_file:
                    for k, v in conflict_list_dict.items():
                        remove_duplicates = list(set(v))
                        if len(remove_duplicates) > 1:
                            conflicting_domain = '\t'.join(remove_duplicates)
                            #writing_conflict_output_file.write(f'{k}\t{conflicting_domain}\n')
                            writing_conflict_output_file.write(f'{k}\n')
                    print(f'Successfully written your {output_conflict_file_name} in directory: {species_folder_in_output}')
                    print('-' * len(output_conflict_file_name) * 5)
                    print('\n')
                                      
###################################################################################################################
        for conflict_list in os.listdir(output_dir):
            if conflict_list.endswith('_conflict.list'):
                species_name = conflict_list.strip().split('_conflict.list')[0]
                conflict_file_path = os.path.join(output_dir, conflict_list)
                output_conflict_bed_name = f'{species_name}_conflict_domain.bed'
                with open(conflict_file_path, 'r') as result_file, open(os.path.join(output_dir, output_conflict_bed_name), 'w') as writing_conflict_bed_file:
                    test_dict = {}
                    for lines in result_file:
                        gene_id_in_conflict_file = lines.strip().split('\t')[0]
                        test_dict[gene_id_in_conflict_file] = []
                        for k, v in conflict_bed_dict.items():
                            #print(f'{k}\t{v}\n')
                            if k == gene_id_in_conflict_file:
                                writing_conflict_bed_file.write(f'{v}')
###################################################################################################################
                                
        for k, v in table_count_dict.items():
            output_table_file_name = f'{k}_counts.txt'
            with open(os.path.join(species_folder_path, output_table_file_name), 'w') as writing_table_file:
                # Write header
                writing_table_file.write(k + '\t'.join([''] + list(v.keys())) + '\n')
        
                # Write species name and counts
                writing_table_file.write(f'{k}\t')
                for count in v.values():
                    writing_table_file.write(f'{count}\t')
                writing_table_file.write('\n')


        # Create Fasta files using '.list' against 'clean.fasta'
        for list_file in os.listdir(species_folder_path):
            if list_file.endswith('.list'):
                list_file_path = os.path.join(species_folder_path, list_file)
                name = list_file.replace('.list', '')
                output_file_name = f'{list_file.split(".list")[0]}.fasta'
                output_file_name_path = os.path.join(species_folder_path, output_file_name)

                for clean_files in os.listdir(directory):
                    if clean_files.endswith('_clean.fasta') and clean_files.startswith(name.split('_')[0]):
                        clean_fasta_path = os.path.join(directory, clean_files)
                        print('-' * len(output_file_name) * 5)
                        print(f'Running: seqkit grep -f {list_file_path} {clean_fasta_path} -o {output_file_name_path}')
                        print('\n')
                        running_seqkit = f'seqkit grep -f {list_file_path} {clean_fasta_path} -o {output_file_name_path}'
                        subprocess.run(running_seqkit, shell = True, check = True)
                print(f'Successfully ran seqkit for {output_file_name} in directory: {species_folder_path}')
                print('-' * len(output_file_name) * 5)
                print('\n')
               
            elif list_file.endswith('_domain.bed'):
                list_file_path = os.path.join(species_folder_path, list_file)
                name = list_file.replace('_domain.bed', '')
                output_file_name = f'{list_file.split("_domain.bed")[0]}_domain.fasta'
                output_file_name_path = os.path.join(species_folder_path, output_file_name)  

                for clean_files in os.listdir(directory):
                    if clean_files.endswith('_clean.fasta') and clean_files.startswith(name.split('_')[0]):
                        clean_fasta_path = os.path.join(directory, clean_files)
                        print('-' * len(output_file_name) * 4)
                        print(f'Running: seqkit subseq --update-faidx --bed {list_file_path} {clean_fasta_path} -o {output_file_name_path}')
                        print('\n')
                        running_seqkit_subseq = f'seqkit subseq --update-faidx --bed {list_file_path} {clean_fasta_path} -o {output_file_name_path}'
                        subprocess.run(running_seqkit_subseq, shell = True, check = True)
                print(f'Successfully ran seqkit for {output_file_name} in directory: {species_folder_path}')
                print('-' * len(output_file_name) * 4)
                print('\n')

        for fasta_file in os.listdir(species_folder_path):
            if fasta_file.endswith('_domain.fasta') and not fasta_file.endswith('_conflict_domain.fasta'):
                fasta_file_path = os.path.join(species_folder_path, fasta_file)
                fasta_name = fasta_file.replace('_domain.fasta', '_muscled_domain.fasta')
                output_muscled_file_name = os.path.join(species_folder_path, fasta_name)
                print('-' * len(fasta_name) * 3)
                print(f'Running: muscle -in {fasta_file_path} -out {output_muscled_file_name}')
                running_muscle = f'muscle -in {fasta_file_path} -out {output_muscled_file_name}'
                subprocess.run(running_muscle, shell = True, check = True)
                print(f'Successfully ran Muscle for {fasta_name} in directory: {species_folder_path}')
                print('-' * len(fasta_name) * 3)
                print('\n')   

        for muscled_file in os.listdir(species_folder_path):
            if muscled_file.endswith('_muscled_domain.fasta'):
                muscled_file_path = os.path.join(species_folder_path, muscled_file)
                trimal_name = muscled_file.replace('_muscled_domain.fasta', '_trimal_muscled_domain.fasta')
                output_trimal_file_name = os.path.join(species_folder_path, trimal_name)
                print('-' * len(trimal_name) * 2)
                print(f'Runnning: trimal -in {muscled_file_path} -out {output_trimal_file_name} -gt 0.50 -cons 60')
                running_trimal = f'trimal -in {muscled_file_path} -out {output_trimal_file_name} -gt 0.50 -cons 60'
                subprocess.run(running_trimal, shell = True, check = True)
                print(f'Successfully ran Trimal for {trimal_name} in directory: {species_folder_path}')
                print('-' * len(trimal_name) * 2)
                print('\n')

        for seqkit_fasta_files in os.listdir(species_folder_path):
            if seqkit_fasta_files.endswith('.fasta'):
                seqkit_file_path = os.path.join(species_folder_path, seqkit_fasta_files)
                running_sed = f"sed -i 's/:.//g' {seqkit_file_path}"
                subprocess.run(running_sed, shell = True, check = True)

    #species_dict = {}
    #for species_file in os.listdir(output_dir):
    #    if species_file.endswith('_domain.fasta')

        if weblogo:
            for muscle_or_trimal_file in os.listdir(species_folder_path):
                if muscle_or_trimal_file.endswith('_muscled_domain.fasta'):
                    muscle_path = os.path.join(species_folder_path, muscle_or_trimal_file)
                    weblogo_muscle_name = muscle_or_trimal_file.replace('_muscled_domain.fasta', '_muscled_domain.pdf')
                    output_weblogo_muscle = os.path.join(species_folder_path, weblogo_muscle_name)
                    print('-' * len(weblogo_muscle_name) * 2)
                    print(f'Running: weblogo -f {muscle_path} -D fasta -o {output_weblogo_muscle} -F pdf --resolution 400')
                    running_weblogo_for_muscle = f'weblogo -f {muscle_path} -D fasta -o {output_weblogo_muscle} -F pdf --resolution 400'
                    subprocess.run(running_weblogo_for_muscle, shell = True, check = True)
                    print(f'Successfully ran Weblogo for {weblogo_muscle_name} in directory: {species_folder_path}')
                    print('-' * len(weblogo_muscle_name) * 2)
                    print('\n')

                elif muscle_or_trimal_file.endswith('_trimal_muscled_domain.fasta'):
                    trimal_path = os.path.join(species_folder_path, muscle_or_trimal_file)
                    weblogo_trimal_name = muscle_or_trimal_file.replace('_trimal_muscled_domain.fasta', '_trimal_muscled_domain.pdf')
                    output_weblogo_trimal = os.path.join(species_folder_path, weblogo_trimal_name)
                    print('-' * len(weblogo_trimal_name) * 2)
                    print(f'Running: weblogo -f {trimal_path} -D fasta -o {output_weblogo_trimal} -F pdf --resolution 400')
                    running_weblogo_for_trimal = f'weblogo -f {trimal_path} -D fasta -o {output_weblogo_trimal} -F pdf --resolution 400'
                    subprocess.run(running_weblogo_for_trimal, shell = True, check = True)
                    print(f'Successfully ran Weblogo for {weblogo_trimal_name} in directory: {species_folder_path}')
                    print('-' * len(weblogo_trimal_name) * 2)
                    print('\n')

def main():
    # Initialize the argument parser
    parser = argparse.ArgumentParser(description="List files from specified directories and set number of CPUs.")
    # Add arguments for input folder, database folder, and number of CPUs
    parser.add_argument("--input", help="Path to the input folder", required=True)
    parser.add_argument("--db", help="Path to the database folder", required=True)
    parser.add_argument("--CPU", help="Number of CPUs", default=2, type=int, required=False)
    parser.add_argument("--logging", help="Enable logging", action='store_true', required=False)
    parser.add_argument("--weblogo", help="Draw weblogo using muscled/aligned sequences", action='store_true')
    
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
        run_hmm(args.input, args.CPU, database_path, args.weblogo)
    else:
        logging.error("Database creation failed or was skipped.")

if __name__ == "__main__":
    main()
