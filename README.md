## GPCR identification

### Environment creation
```bash
conda create -y -n GPCR -c bioconda -c conda-forge -c predector  python=3.11 seqkit hmmer muscle=3.8.1551 weblogo transdecoder signalp6
conda activate GPCR
git clone git@github.com:plantgenomicslab/HMMerMe.git
cd HMMerMe
```
### Run example
usage: run.py [-h] --input INPUT_FOLDER --db DB_FOLDER --CPU CPU [--logging]
```
 python run.py --input Input --db Database --CPU 4
```
output is `output` folder

### To DO Tasks

Within the `output` directory:
1. The `output` directory contains folders named after species. Within each species' folder, there is a file named `{Species_name}.hmm_results`.
2. Omit the first line of `{Species_name}.hmm_results` if it starts with `#`.
3. The fourth column, which lists queries like 7tm_1 and 7tm_2, will be expanded in the future.
4. Generate the following outputs:
   - A file named `{species}_{domain_name=4th_column}.list`, containing only IDs from the first column. Create a separate file for each domain name in the fourth column.
   - A file named `{species}_{domain_name}.fasta`, generated by executing `seqkit grep -f {species}_{domain_name=4th_column}.list INPUT_FOLDER/{species_name}_clean.fasta`.
   - A file named `{species}_{domain_name}_domain.bed`, including data from the first, 18th, and 19th columns. Subtract 1 from the number in the 18th column and append the domain name to the end column.
   - A file named `{species}_{domain_name}_domain.fasta`, created by running `seqkit subseq --update-faidx --bed {species}_{domain_name}_domain.bed INPUT_FOLDER/{species_name}_clean.fasta`.

5. Identify conflicts:
   - Generate `{species}_conflict.list` with Gene ID, 7tm_1, 7tm_2.
   - Create `{species}_conflict.fasta`.
   - Produce `{species}_conflict_domain.bed`, including data from the first, 18th, and 19th columns, append the domain name to the last column, and subtract 1 from the 18th column number.
   - Form `{species}_conflict_domain.fasta` by executing `seqkit subseq --update-faidx --bed {species}_conflict_domain.bed INPUT_FOLDER/{species_name}_clean.fasta`.

6. Construct a table:
   - Use species names as row identifiers and domain names as column headers.
   - Name the file `{OUTPUT}_counts.txt`.
   - the table is under `output` folder.

7. Generate WebLogo:
   - Utilize results from step 4, specifically `{species}_{domain_name}_domain.fasta`.
   - Align sequences using MUSCLE and then create the WebLogo for each species.
   - the output will be under species folder.

8. Combine domain sequences:
   - Combine `{species}_{domain_name}_domain.fasta` files from step 4 into a single file across species.
   - the table is under `output` folder.

9. Script options:
   - Include an `--output` option for specifying the output directory, with `output` as the default.
   - Add a `--CPU` option with a default value of 2.
   - Add an option to enable or disable WebLogo generation, disabled by default.
   - 
