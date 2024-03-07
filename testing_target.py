import os

database_alignment_directory = '/data/gpfs/assoc/pgl/joel/data/HMMerMe/Database_alignment'
directory_path = os.path.join(database_alignment_directory)

output_directory = '/data/gpfs/assoc/pgl/joel/data/HMMerMe/output'
output_path = os.path.join(output_directory)

afa_files_dict = {} #k = Afa file name : v = path to Afa file
for afa_files in os.listdir(directory_path):
    afa_file_name = afa_files.split('.afa')
    afa_file_key = afa_file_name[0]
    afa_files_dict[afa_file_key] = os.path.join(directory_path, afa_files)

for k, v in afa_files_dict.items():
    print(f'{k}\t{v}')
    # GPB2    /data/gpfs/assoc/pgl/joel/data/HMMerMe/Database_alignment/GPB2.afa

for species in os.listdir(output_path):
    print('-----------')
    print(species)
    print('-----------')
    accessing_species_files = os.path.join(output_directory, species)
    for fasta_files in os.listdir(accessing_species_files):
        if fasta_files.endswith('_muscled_domain.fa'):
            for k, v in afa_files_dict.items():
                if k in fasta_files:
                    print(f'muscle -profile -in1 {v} -in2 {fasta_files} -out {fasta_files.replace("_muscled_domain.fasta", "_muscled_combined_domain.afa")}')