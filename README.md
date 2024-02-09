## GPCR identification

### Environment creation
```bash
conda create -y -n GPCR -c bioconda -c conda-forge python=3.11 seqkit hmmer muscle=3.8.1551 weblogo
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

### To DO
From `output` 
1. the output folder includes species name folder, in the species name folder it will have a file {Species_name}.hmm_results
2. From {Species_name}.hmm_results, skip the first line start with #
3. 4th column is query such as 7tm_1 and 7tm_2. It will be added more later.
4. Make the output as below
- {species}_{domain_name=4th_column}.list. This will include 1st column ID only. Each 4th column need to have a individual file.
- {species}_{domain_name}.fasta. seqkit grep -f  {species}_{domain_name=4th_column}.list INPUT_FOLDER/{species_name}_clean.fasta
- {species}_{domain_name}_domain.bed. This will include 1st, 18th and 19th coulumn. Then 18th column need to (number substract 1), add domain name to the end column 
- {species}_{domain_name}_domain.fasta seqkit subseq --update-faidx --bed {species}_{domain_name}_domain.bed INPUT_FOLDER/{species_name}_clean.fasta

5. Find the conflict
- {species}_conflict.list -> Gene ID, 7tm_1, 7tm_2
- {species}_conflict.fasta 
- {species}_conflict_domain.bed This will include 1st, 18th and 19th coulumn. add domain name to the end column. Then 18th column need to (number substract 1)
- {species}_conflict_domain.fasta -> seqkit subseq --update-faidx --bed {species}_conflict_domain.bed INPUT_FOLDER/{species_name}_clean.fasta

6. make a table
- Species name is row name, domain name column name
- file name is {OUTPUT}_counts.txt

7. Weblogo
- use results from 4 :{species}_{domain_name}_domain.fasta
- Align with muscle then draw the weblogo  

8. Merge domain sequences
 - merge {species}_{domain_name}_domain.fasta from 4, merge them into one across species.

9. add option for script
- --output option for script default is output
- add --CPU default is 2
- add option for Weblogo or not default is not


