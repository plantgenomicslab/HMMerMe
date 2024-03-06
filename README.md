## Installation
### Create 'HMM analysis' environment with the follwoing software tools installed

`conda create -y -n HMM -c bioconda -c conda-forge python=3.11 seqkit hmmer muscle=3.8.1551 weblogo transdecoder easel diamond trimal pymsaviz`

Activating said environment

`conda activate HMM`

### Conda installation
Please refer below link
https://docs.anaconda.com/free/miniconda/[https://docs.anaconda.com/free/miniconda/]

### Download and utlize 'HMMer'

In your working directory of choice

`git clone git@github.com:plantgenomicslab/HMMerMe.git`

If you have not created an SSH key, please do so:

Go into 'HMM' folder

`cd HMMerME`


## Running 'HMMerME'

Main code will be in `run.py`

`run.py` takes in 2 required arguments: `--input` and `--db`

`--input` will target the `Input` directory that has your `.fasta` or `.fas` files

`--db` will target the 'Database' directory that has your '.hmm' files, called after an HMM search

`--CPU` is one of the options that utilizes a specific core. If not called, default is set to use 2 cores (4 threads)

`--visualization` is another option that calls 'weblogo' and 'pymsaviz' to visualize your data

To run the command, go into the 'HMMerME' directory as the current working directory as simply type this:

`python run.py --input Input/ --db Database/ [OPTIONS --CPU, --visualization]`

Remember, if you want a weblogo and visualized aligned sequences, run `--visualization`

`python run.py --input Input/ --db Database/ --visualization`

### Input and Output Files

Input files include

```
Database/{Homology_domain}.hmm
Input/{Species}.fasta
```

Output files include

```
{Species}.hmm_results
{Species}_{Homology_domain}.list
{Species}_{Homology_domain}.fasta
{Species}_{Homology_domain}_domain.bed
{Species}_{Homology_domain}_domain.fasta

{Species}_conflict.list
{Species}_conflict.fasta
{Species}_conflict_domain.bed
{Species}_conflict_domain.fasta

{Species}_counts.txt

{Species}_{Homology_domain}_muslced_domain.fa
{Species}_{Homology_domain}_muscled_trimal_domain.fa

{Species}_{Homology_domain}_muscled_combined_domain.afa
{Species}_{Homology_domain}_muscled_combined_trimal_domain.afa

{Species}_{Homology_domain}_muscled_combined_trimal_domain.pdf
{Species}_{Homology_domain}_muscled_combined_trimal_domain.png
```
## Quick explanation

Within the `Database` directory, there are `.hmm` files. These files are your HMM profile files. An HMM profile is simply a statistical algorithms that detects homology between your query and a predefined large sequence databases (Johnson et al. 2010). These HMM profile files are found in the `Database` directory. Additonally, within your `Input` directory, there are your species fasta files. Having both the fasta files and HMM profile files will allow `run.py` to find the homology sequence between your fasta files and HMM profiles. 

Lets take a look at the format for species AaegyptiLVPWY. There are exactly 28,392 genes within the AaegyptiLVPWY fasta file. In the Database folder, there are 49 HMM profile that were found from calling HMM search. `run.py` will try to find and match all possible homology between your fasta file and the HMM profile. If found, you should obtain this type of format: {Species}_{Homology_domain}. Additionally, there are a bunch of output files called after `run.py`. Let go over them in breif

The first file that outputs is the {Species}.hmm_results. After Processing each of your FASTA file against a large datbase called by HMMER, a homology results will show with an E-value of less thatn 1e-10 (Remember, lower the E-value, greater the significance that the search was not due to chance but for biological significance). Within the `.hmm_results` file, Column 1 represents the domain gene, Column 2 represents the domain name, Column 18 and 19 represents the alignment coordinate form and to respectivley. Your `.hmm_search` file is important as it will be used to model your `.list`, `.fasta`, `domain.bed`, and `domain.fasta` files. 

The `.list`, `.fasta`, `domain.bed`, and `domain.fasta` files have the formats: {Species}_{Homology_domain} where Species is your species of interest and the name of the Homologous domain group. `.list` will tel you the domain genes that is represented from the search, `.fasta`, will search domain gene sequences that are found in the the fasta files, `domain.bed` will give you a bed format that includes where the domain gene is found in coordinate format, and `domain.fasta` will print out respective subsets of domain gene sequences in fasta format form by using `domain.bed`.

From creating the `.list`, `.fasta`, `domain.bed`, and `domain.fasta` files, it will also create `conflict.list`, `conflict.fasta`, `conflict_domain.bed`, and `conflict_domain.fasta` files. To put it in simple terms, some domain genes can have the same sequences but differnt domain names. Therefore, the `conflict.list` file will represents the domain genes that have conflicing domain names, `conflict.fasta` will output the fasta sequences of that conflicting genes, `conflict_domain.bed` will give you a bed format of the conflicting domain gene in coordinate format, and `conflict.fasta` will extract subsets of conflciting gene sequences using `conflict_domain.bed`

Next is the `counts.txt`. This file will give you the total amount of Domain genes that were distinguished from your `.hmm_search` file.
| AaegyptiLVPWY | RYamideLuqin | Prothoracicotropichormone | SIfamide | CCHamide1 |
|---------------|--------------|---------------------------|----------|-----------|
| AaegyptiLVPWY | 1            | 3                         | 1        | 1         |

The `muslced_domain.fa` files are simply aligned fasta sequences using the Muscle program version 3.8.1551 (Edgar et al 2021). Additionally, some sequences have to much of a gap, we will call trimal to remove these gaps in a file format of `muscled_trimal_domain.fa`.

The `muscled_combined_domain.afa`files combines the sequence alignments of the HMM search files in the form of `.afa` files (.fasta format). In turn, `run.py` will try to combine your `muscled_combined_domain.fa` and the `{Homology_domain}.afa` files in the `Database_alignment` directory if a match between the Homology name is found. the 'muscled_combined_trimal_domain.afa' simply calls trimal to remove unecessary gaps.

Finally `.pdf` and `.png` files are made. If you called `--visualziation` it will call Weblogo v3 (Crooks et al 2004) and Pymsaviz (https://pypi.org/project/pyMSAviz/) to create the `.pdf` and `.png` files respectively. Here are some visualization:

[AaegyptiLVPWY_Adipokinetic_Corazonin_muscled_combined_trimal_domain.pdf](https://github.com/plantgenomicslab/HMMerMe/files/14514243/AaegyptiLVPWY_Adipokinetic_Corazonin_muscled_combined_trimal_domain.pdf)

![AaegyptiLVPWY_Adipokinetic_Corazonin_muscled_combined_trimal_domain](https://github.com/plantgenomicslab/HMMerMe/assets/137996393/75e025cf-6219-40db-9465-78a11d80f7c4)


   - Include an `--output` option for specifying the output directory, with `output` as the default.
   - Add a `--CPU` option with a default value of 2.
   - Add an option to enable or disable WebLogo generation, disabled by default.

## Custom database creation
muscle -in {Homology_domain}.fasta  -out {Homology_domain}.aln -clw
trimal -in {Homology_domain}.aln -out {Homology_domain}_trimmed.aln -gt 0.50 -cons 60
esl-reformat stockholm  {Homology_domain}_trimmed.aln > {Homology_domain}.sto
hmmbuild  {Homology_domain}.hmm {Homology_domain}.sto

