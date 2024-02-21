import argparse, subprocess, re, os 

parser = argparse.ArgumentParser(description = 'File test')
parser.add_argument('-list', '--list', action = 'store_true', help = 'Option: Create a "List" File. Contains: ID\tDomain', required = False)
parser.add_argument('-fasta', '--fasta', action = 'store_true', help = 'Option: Create a "Fasta" File. Uses Seqkit', required = False)
parser.add_argument('-bed', '--bed', action = 'store_true', help = 'Option: Create a "BED" File. Contains: ID\tfrom\tto', required = False)
parser.add_argument('-subseq', '--subseq', action = 'store_true', help = 'Option: Create a "Fasta" FIle using a clean fasta file.', required = False)
parser.add_argument('-conflict_list', '--conflict_list', action = 'store_true', help = 'Option: Create a "Conflict" File list.', required = False)
parser.add_argument('-conflict_id', '--conflict_id', action = 'store_true', help = 'Option: Create a ID file for seqkit command', required = False)
parser.add_argument('-conflict_fasta', '--conflict_fasta', action = 'store_true', help = 'Option: Create a "Conflict" fasta file', required = False)
parser.add_argument('-conflict_bed', '--conflict_bed', action = 'store_true', help = 'Option: Create a "Conflict" BED file', required = False)
parser.add_argument('-conflict_clean', '--conflict_clean', action = 'store_true', help = 'Option: Create a "Clean" fasta conflict file', required = False)
parser.add_argument('-table', '--table', action = 'store_true', help = 'Option: Create a "Table" count File', required = False)
parser.add_argument('-muscle', '--muscle', action = 'store_true', help = 'Option: Use Muscle to align sequences', required = False)
parser.add_argument('-trimal', '--trimal', action = 'store_true', help = 'Option: Pherhaps your WebLogo design is not correct, introcuding trimal!', required = False)
parser.add_argument('-weblogo', '--weblogo', action = 'store_true', help = 'Option: Create a WebLogo design', required = False)
parser.add_argument('-fix', '--fix', action = 'store_true', help = 'Option: For files that used subseq, fix the format', required = False)
parser.add_argument('-combine_domain', '--combine_domain', action = 'store_true', help = 'Option: Combine the domain sequences for species', required = False)
parser.add_argument('-i', '--input_file', type = str, help = 'Input file', required = False)
parser.add_argument('-c', '--clean_file', type = str, help = 'Clean file', required = False)
parser.add_argument('-o', '--output_file', type = str, help = 'Output File name', required = False)
parser.add_argument('-d', '--output_dir', type = str, help = 'Directory of Output File destination', required = False)
args = parser.parse_args()

# Functions for Step 4
def domain_list(input_file, output_file, output_dir):
    if args.list:
        list_file_dict = {}
        with open(input_file, 'r') as file:
            for lines in file:
                if lines.startswith('#') or lines.startswith('-'):
                    continue
                else:
                    tabbed_lines = re.sub(r'\s+', '\t', lines)
                    columns = tabbed_lines.strip().split('\t')
                    id, domain = columns[0], columns[3]
                    if domain not in list_file_dict:
                        list_file_dict[domain] = []
                    else:
                        list_file_dict[domain].append(f'{id}')

        for k, v in list_file_dict.items():
            if output_dir:
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                    if not output_file:
                        output_file = f'{os.getcwd()}_{k}.list'
                    output_file_path = os.path.join(output_file)
            else:
                output_file_path = output_file or f'{os.getcwd()}_{k}.list'

            with open(output_file_path, 'w') as output_file_destination:
                for values in v:
                    output_file_destination.write(f'{values}\n')
                print('\n')
                print(f'--------------------------------------------------------------------------------------------\n')
                print(f'Success! Written your Domain list to:\n')
                print(f'{output_file_destination}\n')
                print(f'--------------------------------------------------------------------------------------------\n')
                print('\n')

def fasta(input_file, clean_file):
    if args.fasta:
        input_id_file = os.path.join(os.getcwd(), input_file)
        input_clean_file = os.path.join(os.getcwd(), clean_file)
        output_file_name = input_id_file.replace(".list", ".fasta")
        subprocess.run(['seqkit', 'grep', '-f', f'{input_id_file}', f'{input_clean_file}', '-o', f'{output_file_name}'])
        print('\n')
        print(f'------------------------------------------------------------------------------------\n')
        print(f'Success! Written your FASTA file to:\n')
        print(f'{os.path.join(os.getcwd(), output_file_name)}\n')
        print(f'------------------------------------------------------------------------------------\n')
        print('\n')

def bed(input_file, output_file, output_dir):
    if args.bed:
        bed_file_dict = {}
        with open(input_file, 'r') as file:
            for lines in file:
                if lines.startswith('#') or lines.startswith('-'):
                    continue
                else:
                    tabbed_lines = re.sub(r'\s+', '\t', lines)
                    columns = tabbed_lines.strip().split('\t')
                    id, domain, ali_from, ali_to = columns[0], columns[3], columns[17], columns[18]
                    if domain not in bed_file_dict:
                        bed_file_dict[domain] = []
                    else:
                        bed_file_dict[domain].append(f'{id}\t{int(ali_from) - 1}\t{ali_to}')

        for k, v in bed_file_dict.items():
            if output_dir:
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                    if not output_file:
                        output_file = f'{os.getcwd()}_{k}.bed'
                    output_file_path = os.path.join(output_file)
            else:
                output_file_path = output_file or f'{os.getcwd()}_{k}.bed'

            with open(output_file_path, 'w') as output_file_destination:
                for values in v:
                    output_file_destination.write(f'{values}\n')
                    print('\n')
                    print(f'-----------------------------------------------------------------------------------------\n')
                    print(f'Success! Written your BED file to:\n')
                    print(f'{output_file_destination}\n')
                    print(f'-----------------------------------------------------------------------------------------\n')
                    print('\n')

def clean_fasta(input_file, clean_file):
    if args.subseq:
        input_id_file = os.path.join(os.getcwd(), input_file)
        input_clean_file = os.path.join(os.getcwd(), clean_file)
        output_file_name = input_id_file.replace(".bed", "_domain.fasta")
        subprocess.run(['seqkit', 'subseq', '-R', '-U', '--bed', f'{input_id_file}', f'{input_clean_file}', '-o', f'{output_file_name}'])
        print('\n')
        print(f'-------------------------------------------------------------------------------------------------------------\n')
        print(f'Success! written your Updated FASTA sequence using "Subseq" to:\n')
        print(f'{os.path.join(os.getcwd(), output_file_name)}\n')
        print(f'-------------------------------------------------------------------------------------------------------------\n')
        print('\n')

# Functions for Step 5
def create_conflict_list(input_file, output_file, output_dir):
    if args.conflict_list:
        if output_dir:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
                if not output_file:
                    output_file = f'{os.getcwd()}_conflict.list'
                output_file_path = os.path.join(output_file)
        else:
            output_file_path = output_file or f'{os.getcwd()}_conflict.list'

        list_file_dict = {}
        with open(input_file, 'r') as file:
            for lines in file:
                if lines.startswith('#') or lines.startswith('-'):
                    continue
                else:
                    tabbed_lines = re.sub(r'\s+', '\t', lines)
                    columns = tabbed_lines.strip().split('\t')
                    id, domain = columns[0], columns[3]
                    if id not in list_file_dict:
                        list_file_dict[id] = []
                    list_file_dict[id].append(domain)
        
        #for k, v in list_file_dict.items():
            #print(f'{k}\t{v}\n')

        with open(output_file_path, 'w') as output_file_destination:
            for k, v in list_file_dict.items():
                unique_values = list(set(v))
                if len(unique_values) > 1:
                    conflicting_domain = "\t".join(unique_values)
                    output_file_destination.write(f'{k}\t{conflicting_domain}\n')
                    print('\n')
                    print(f'---------------------------------------------------------------------------------------------------\n')
                    print(f'Success! Written your Conflict list file to:\n')
                    print(f'{output_file_destination}\n')
                    print(f'---------------------------------------------------------------------------------------------------\n')
                    print('\n')

def create_conflict_IDs(input_file, output_file, output_dir):
    if args.conflict_id:
        if output_dir:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
                if not output_file:
                    output_file = f'{os.getcwd()}_ID_conflict.fasta'
                output_file_path = os.path.join(output_file)
        else:
            output_file_path = output_file or f'{os.getcwd()}_ID_conflict.fasta'

        list_file_dict = {}
        with open(input_file, 'r') as file:
            for lines in file:
                if lines.startswith('#') or lines.startswith('-'):
                    continue
                else:
                    tabbed_lines = re.sub(r'\s+', '\t', lines)
                    columns = tabbed_lines.strip().split('\t')
                    id, domain = columns[0], columns[3]
                    if id not in list_file_dict:
                        list_file_dict[id] = []
                    list_file_dict[id].append(domain)

        with open(output_file_path, 'w') as output_file_destination:
            for k, v in list_file_dict.items():
                unique_values = list(set(v))
                if len(unique_values) > 1:
                    output_file_destination.write(f'{k}\n')
                    print('\n')
                    print(f'---------------------------------------------------------------------------------------------------\n')
                    print(f'Success! Written your Conflict list file to:\n')
                    print(f'{output_file_destination}\n')
                    print(f'---------------------------------------------------------------------------------------------------\n')
                    print('\n')

def create_conflict_fasta_list(input_file, clean_file):
    if args.conflict_fasta:
        input_id_file = os.path.join(os.getcwd(), input_file)
        input_clean_file = os.path.join(os.getcwd(), clean_file)
        output_file_name = input_id_file.replace("_ID_conflict.fasta", "_conflict.fasta")
        subprocess.run(['seqkit', 'grep', '-f', f'{input_id_file}', f'{input_clean_file}', '-o', f'{output_file_name}'])
        print('\n')
        print(f'---------------------------------------------------------------------------------------------\n')
        print(f'Success! Written your Conflict FASTA file to:\n')
        print(f'{os.path.join(os.getcwd(), output_file_name)}\n')
        print(f'---------------------------------------------------------------------------------------------\n')
        print('\n')

def create_conflict_bed(input_file, output_file, output_dir):
    if args.conflict_bed:
        bed_file_dict = {}

        if output_dir:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
                if not output_file:
                    output_file = f'{os.getcwd()}_conflict.bed'
                output_file_path = os.path.join(output_file)
        else:
            output_file_path = output_file or f'{os.getcwd()}_conflict.bed'

        with open(input_file, 'r') as file:
            for lines in file:
                if lines.startswith('#') or lines.startswith('-'):
                    continue
                else:
                    tabbed_lines = re.sub(r'\s+', '\t', lines)
                    columns = tabbed_lines.strip().split('\t')
                    id, domain, ali_from, ali_to = columns[0], columns[3], columns[17], columns[18]
                    if domain not in bed_file_dict:
                        bed_file_dict[domain] = []
                    else:
                        bed_file_dict[domain].append(f'{id}\t{int(ali_from) - 1}\t{ali_to}')

        with open(output_file_path, 'w') as output_file_destination:
            for k, v in bed_file_dict.items():
                print(f'{k}\t{v}\n')
                for values in v:
                    output_file_destination.write(f'{values}\t{k}\n')
                print('\n')
                print(f'--------------------------------------------------------------------------------------------------\n')
                print(f'Success! Written your Conflict BED file to:\n')
                print(f'{os.path.join(os.getcwd(), output_file_destination)}\n')
                print(f'--------------------------------------------------------------------------------------------------\n')
                print('\n')

def create_clean_conflict_fasta(input_file, clean_file):
    if args.conflict_clean:
        input_id_file = os.path.join(os.getcwd(), input_file)
        input_clean_file = os.path.join(os.getcwd(), clean_file)
        output_file_name = input_id_file.replace(".bed", "_domain.fasta")
        subprocess.run(['seqkit', 'subseq', '--update-faidx', '--bed', f'{input_id_file}', f'{input_clean_file}', '-o', f'{output_file_name}'])
        print('\n')
        print(f'--------------------------------------------------------------------------------------------------------------------\n')
        print(f'Success! Written your updated Conflict FASTA file using "Subseq" to:\n')
        print(f'{os.path.join(os.getcwd(), output_file_name)}\n')
        print(f'--------------------------------------------------------------------------------------------------------------------\n')
        print('\n')

# Function for Step 6
def table(output_dir, output_file):
    if args.table:

        if output_dir:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
                if not output_file:
                    output_file = f'table.txt'
                output_file_path = os.path.join(output_file)
        else:
            output_file_path = output_file or f'table.txt'

        target_directory = os.getcwd()
        path_to_list_files = os.listdir(target_directory)

        species_count_dict = {}
        for files in path_to_list_files:
            if files.endswith('.hmm_results'):
                with open(os.path.join(target_directory, files), 'r') as result_file:
                    for lines in result_file:
                        if lines.startswith('#') or lines.startswith('-'):
                            continue
                        else:
                            create_tabbed_lines = re.sub(r'\s+', '\t', lines)
                            columns = create_tabbed_lines.strip().split('\t')
                            id, domain = columns[0], columns[3]
                            speciesname = os.path.splitext(files)[0]
                            if speciesname not in species_count_dict:
                                species_count_dict[speciesname] = {}
                            if domain not in species_count_dict[speciesname]:
                                species_count_dict[speciesname][domain] = 1
                            else:
                                species_count_dict[speciesname][domain] += 1
        
        #for key, value in species_count_dict.items():
            #print(f'{key}\t{value}\n')

        all_domains = set()
        for domain_name in species_count_dict.values():
            all_domains.update(domain_name)
        print(all_domains)

        sorted_domain = sorted(all_domains)
        with open(os.path.join(output_file_path), 'w') as output:
            output.write('Species\t')
            for names_of_domain in sorted_domain:
                output.write(f'{names_of_domain}\t')
            output.write('\n')
            for k, v in species_count_dict.items():
                output.write(f'{k}\t')
                for counts in all_domains:
                    total = v.get(counts, 0)
                    output.write(f'{total}\t')
                output.write('\n')
            print('\n')
            print(f'------------------------------------------------------------------------------\n')
            print(f'Success! Formed a Domain Count Table to:\n')
            print(f'{os.getcwd()}\n')
            print('as:\n')
            print(f'table.txt\n')
            print(f'------------------------------------------------------------------------------\n')
            print('\n')

# Function for muscle
def muscle():
    if args.muscle:
        fasta_folder = os.getcwd()
        fasta_files = []
        for file in os.listdir(fasta_folder):
            if file.endswith('_domain.fasta') and not file.endswith('_conflict.fasta'):
                fasta_files.append(file)

        for file in fasta_files:
            input_file_path = os.path.join(fasta_folder, file)
            output_file_path = os.path.join(fasta_folder, f'muscled_{file}')
            subprocess.run(['muscle', '-in', f'{input_file_path}', '-out', f'{output_file_path}'])
            print('\n')
            print(f'-------------------------------------------------------------------------------------------\n')
            print(f'Success! Aligned your {input_file_path} to:\n')
            print(f'{os.getcwd()}\n')
            print('as:\n')
            print(f'muscled_{file}')
            print(f'-------------------------------------------------------------------------------------------\n')
            print('\n')
        print('Thanks for using Muscle!')
        print('\n')

# Function for trimal
def trimal():
    if args.trimal:
        current_fasta_folder = os.getcwd()
        current_files = []
        for file in os.listdir(current_fasta_folder):
            if file.endswith('_domain.fasta') and file.startswith('muscled_') and not file.endswith('_conflict.fasta'):
                current_files.append(file)

        for file in current_files:
            input_file_path = os.path.join(current_fasta_folder, file)
            output_file_path = os.path.join(current_fasta_folder, f"{file.replace('_domain.fasta', '_domain_clean.fasta')}")
            subprocess.run(['trimal', '-in', f'{input_file_path}', '-out', f'{output_file_path}', '-gappyout'])
            print('\n')
            print(f'-------------------------------------------------------------------------------------------\n')
            print(f'Success! Removed unwanted gaps of {input_file_path}. Saved file to:\n')
            print(f'{os.getcwd()}\n')
            print('as:\n')
            print(f"{file.replace('_domain.fasta', '_domain_clean.fasta')}")
            print(f'-------------------------------------------------------------------------------------------\n')
            print('\n')
        print('Thanks for using Muscle!')
        print('\n')

# Function for Step 7
def weblogo():
    if args.weblogo:
        muscled_folder = os.getcwd()
        muscled_files = []
        for file in os.listdir(muscled_folder):
            if file.startswith('muscled_') and file.endswith('_domain.fasta') or file.endswith('_domain_clean.fasta'):
                muscled_files.append(file)
        
        for file in muscled_files:
            input_file_path = os.path.join(muscled_folder, file)
            output_file_path = os.path.join(muscled_folder, f'weblogo_{file.replace(".fasta", ".pdf")}')
            subprocess.run(['weblogo', '-f', f'{input_file_path}', '-D', 'fasta', '-o', f'{output_file_path}', '-F', 'pdf', '--resolution', '400'])
            print('\n')
            print(f'-------------------------------------------------------------------------------------------------\n')
            print(f'Success! Created your Weblogo Design as a PDf to:\n')
            print(f'{os.getcwd()}\n')
            print(f'as:\n')
            print(f'weblogo_{file.replace(".fasta", ".pdf")}\n')
            print(f'-------------------------------------------------------------------------------------------------\n')
            print('\n')

# Function format "_domain.fasta" files correctly
def fix_format(output_file, output_dir):
    if args.fix:
        domain_fasta = []
        domain_file_folder = os.getcwd()
        for file in os.listdir(domain_file_folder):
            if file.endswith('_domain.fasta') and not file.endswith('_conflict_domain.fasta') or file.startswith('muscled_'):
                domain_fasta.append(file)

        for domainFiles in domain_fasta:
            domain_file_path = os.path.join(domain_file_folder, domainFiles)
            remove_species_name = re.sub(r'^[^_]*_', '', domain_file_path)
            domain_name = remove_species_name.replace('_domain.fasta', '')
            if output_dir:
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                    if not output_file:
                        output_file = f'fixed_{domainFiles}'
                    output_file_path = os.path.join(output_file)
            else:
                output_file_path = output_file or f'fixed_{domainFiles}'
            
            with open(domain_file_path, 'r') as domain_files, open(output_file_path, 'w') as output:
                for seqkit_headers_lines in domain_files:
                    if seqkit_headers_lines.startswith('>'):
                        new_header = re.sub(r':.', f'-{domain_name}', seqkit_headers_lines)
                        output.write(new_header)
                    else:
                        output.write(seqkit_headers_lines)

# Function for Step 8
def domain_sequence(output_dir, output_file):
    if args.combine_domain:
        if args.combine_domain:
            species_dict = {}
            species_folder = os.getcwd()
            for species_file in os.listdir(species_folder):
                if species_file.startswith('fixed_') and species_file.endswith('_domain.fasta') and not species_file.endswith('_conflict_domain.fasta') and not species_file.startswith('fixed_muscled_'):
                    species_name = species_file.strip().split('_')[1]

                    if species_name not in species_dict:
                        species_dict[species_name] = []

                    with open(species_file, 'r') as fasta_information:
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
                    output_file_path = os.path.join(output_dir, f'{k}_combined.fasta')
                else:
                    output_file_path = output_file or f'{k}_combined.fasta'

                with open(output_file_path, 'w') as output:
                    for seq in v:
                        output.write(f'{seq}\n')
                    print(f'{output_file_path}')

# Call Functions for Step 4
domain_list(args.input_file, args.output_file, args.output_dir)
fasta(args.input_file, args.clean_file)
bed(args.input_file, args.output_file, args.output_dir)
clean_fasta(args.input_file, args.clean_file)

# Call Functions for Step 5
create_conflict_list(args.input_file, args.output_file, args.output_dir)
create_conflict_IDs(args.input_file, args.output_file, args.output_dir)
create_conflict_fasta_list(args.input_file, args.clean_file)
create_conflict_bed(args.input_file, args.output_file, args.output_dir)
create_clean_conflict_fasta(args.input_file, args.clean_file)

# Call Functions for Step 6
table(args.output_dir, args.output_file)
muscle()
trimal()
weblogo()
fix_format(args.output_dir, args.output_file)
domain_sequence(args.output_dir, args.output_file)