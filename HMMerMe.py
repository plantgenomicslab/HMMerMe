#!/usr/bin/env python

import os
import subprocess
from Bio import SeqIO

# Set umask
os.umask(0o007)

# Define input directories
input_transcriptomes_dir = 'input_Transcriptomes'
input_hmm_profiles_dir = 'input_HMM_profiles'

# Define output directories
output_dir = 'Nailed_it'
conflicted_output_dir = os.path.join(output_dir, 'conflicted_HMMer_output')
processed_input_dir = 'Processed_input_files/Transcriptomes'
database_hmm_dir = 'Database_with_HMM_profiles'

# Create necessary directories
os.makedirs(output_dir, exist_ok=True)
os.makedirs(conflicted_output_dir, exist_ok=True)
os.makedirs(processed_input_dir, exist_ok=True)
os.makedirs(database_hmm_dir, exist_ok=True)

# Function to remove special characters and shorten FASTA headers
def clean_sequence_header(header):
    # Remove special characters from the header
    header = ''.join(char if char.isalnum() or char in ['-', '_', '|'] else ' ' for char in header)
    # Shorten the header by keeping only the part before the first space
    header = header.split()[0]
    return header

# Level 1: Iterate through input transcriptomes
for transcriptome in os.listdir(input_transcriptomes_dir):
    if transcriptome.endswith(('.fasta', '.fas', '.fa')):
        species = os.path.splitext(transcriptome)[0]

        print(f"Starting with {species}")

        transcriptome_path = os.path.join(input_transcriptomes_dir, transcriptome)

        # Read and process the transcriptome using Biopython
        records = []
        for record in SeqIO.parse(transcriptome_path, "fasta"):
            # Clean the sequence header
            record.id = clean_sequence_header(record.id)
            records.append(record)

        # Write the cleaned transcriptome to a temporary file
        temp_transcriptome_path = os.path.join('temp_transcriptomes', transcriptome)
        os.makedirs(os.path.dirname(temp_transcriptome_path), exist_ok=True)
        SeqIO.write(records, temp_transcriptome_path, "fasta")

        # Create temporary directories
        tmp_conflicts_dir = os.path.join('tmp_conflicts_dir', species)
        os.makedirs(tmp_conflicts_dir, exist_ok=True)
        os.makedirs(os.path.join(tmp_conflicts_dir, 'Identifiers_separated'), exist_ok=True)
        os.makedirs(os.path.join(tmp_conflicts_dir, 'Identifiers_combined'), exist_ok=True)
        os.makedirs(os.path.join(tmp_conflicts_dir, 'Sequences'), exist_ok=True)
        os.makedirs(os.path.join(tmp_conflicts_dir, 'tmp_IDs'), exist_ok=True)
        os.makedirs(os.path.join(tmp_conflicts_dir, 'Unique_identifiers'), exist_ok=True)
        os.makedirs(os.path.join(tmp_conflicts_dir, 'Unique_tmp_IDs'), exist_ok=True)

        # Level 2: Iterate through HMM profiles
        for hmm_profile in os.listdir(input_hmm_profiles_dir):
            if hmm_profile.endswith('.hmm'):
                hmm_profile_path = os.path.join(input_hmm_profiles_dir, hmm_profile)
                hmm_profile_name = os.path.splitext(hmm_profile)[0]
                prediction_file = f"{species}__{hmm_profile_name}_predictions.fasta"
                identifier_list = f"{species}__{hmm_profile_name}_identifiers.fasta"
                initial_hmmer_results = f"{species}__{hmm_profile_name}_evalues.fasta"
                output1 = f"{species}__{hmm_profile_name}_output1.fasta"
                output2 = f"{species}__{hmm_profile_name}_output2.fasta"
                output3 = f"{species}__{hmm_profile_name}_output3.fasta"
                output4 = f"{species}__{hmm_profile_name}_output4.fasta"

                print(f"Search for {hmm_profile_name}")

                # Perform hmmsearch
                subprocess.run(['hmmsearch', '-E', '1e-10', hmm_profile_path, transcriptome_path, '>', f'{tmp_conflicts_dir}/{output1}'], shell=True, check=True)

                # Process the output
                with open(f'{tmp_conflicts_dir}/{output1}', 'r') as f:
                    lines = f.readlines()
                
                # Extract relevant lines
                start_index = None
                end_index = None
                for i, line in enumerate(lines):
                    if '--- full sequence ---' in line:
                        start_index = i + 1
                    elif 'Domain annotation for each sequence (and alignments)' in line:
                        end_index = i - 1
                        break
                
                if start_index is not None and end_index is not None:
                    relevant_lines = lines[start_index:end_index + 1]
                    output_content = '\n'.join(relevant_lines)
                    
                    # Save the processed output
                    with open(f'Nailed_it/original_HMMer_output/{initial_hmmer_results}', 'w') as f:
                        f.write(output_content)
                
                # Continue with other processing steps (not included in this code snippet)

        # Move the processed transcriptome to the "processed" folder
        transcriptome_content = transcriptome_content.replace('___', '|').replace('  ', ' ')
        with open(f'Processed_input_files/Transcriptomes/{transcriptome}', 'w') as f:
            f.write(transcriptome_content)

# Move HMM profiles to the "Database_with_HMM_profiles" folder
subprocess.run(['mv', 'input_HMM_profiles/*', 'Database_with_HMM_profiles/'], shell=True, check=True)

# Count and output the number of predicted genes
subprocess.run(['grep', '-c', '>', './Nailed_it/original_HMMer_output/*_predictions.fasta', '>', './Nailed_it/original_HMMer_output/00_Number_of_all_predicted_genes.txt'], shell=True, check=True)
subprocess.run(['sed', '-i', 's,./Nailed_it/original_HMMer_output/,,g; s,__,\t,g; s,_predictions.fasta:,\t,g', './Nailed_it/original_HMMer_output/00_Number_of_all_predicted_genes.txt'], shell=True, check=True)

subprocess.run(['grep', '-c', '>', './Nailed_it/conflicted_HMMer_output/*_conflicts.fasta', '>', './Nailed_it/conflicted_HMMer_output/00_Number_of_conflicted_predicted_genes.txt'], shell=True, check=True)
subprocess.run(['sed', '-i', 's,./Nailed_it/conflicted_HMMer_output/,,g; s,__,\t,g; s,_conflicts.fasta:,\t,g', './Nailed_it/conflicted_HMMer_output/00_Number_of_conflicted_predicted_genes.txt'], shell=True, check=True)

subprocess.run(['grep', '-c', '>', './Nailed_it/*_no_conflicts.fasta', '>', './Nailed_it/00_Number_of_predicted_genes_without_conflicts.txt'], shell=True, check=True)
subprocess.run(['sed', '-i', 's,./Nailed_it/,,g; s,__,\t,g; s,_no_conflicts.fasta:,\t,g', './Nailed_it/00_Number_of_predicted_genes_without_conflicts.txt'], shell=True, check=True)

print("--------------------------------------------------------------")
print(" ")
print("Done!")
print(" ")
print(".hmm input profiles were moved to 'Database_with_HMM_profiles'")
print(" ")
print("Input Transcriptomes were moved to 'Processed_input_files/Transcriptomes'")
print(" ")
print("Original HMMer predictions are in 'Nailed_it/original_HMMer_output'")
print(" ")
print("Sequences that were predicted by more than one HMM profile (= conflicts) are in 'Nailed_it/conflicted_HMMer_output'")
print(" ")
print("Sequences that were only predicted once (= no conflicts) are in 'Nailed_it'")
print(" ")
print("Number of predicted sequences (original, conflicts, no-conflicts) are listed in the 00_Number...txt files in the corresponding folders")
