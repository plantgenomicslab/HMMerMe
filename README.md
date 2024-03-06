## Setup Instructions

### Creating the 'HMM Analysis' Environment

Install the required tools in a new Conda environment named "HMM" using the following command:

```bash
conda create -y -n HMM -c bioconda -c conda-forge python=3.11 seqkit hmmer muscle=3.8.1551 weblogo transdecoder easel diamond trimal pymsaviz
```

To activate the newly created environment, execute:

```bash
conda activate HMM
```

### Installing Conda

For Conda installation, visit the following [link](https://docs.anaconda.com/free/miniconda/).

### Setting Up 'HMMer'

Select a desired working directory and clone the HMMerMe repository:

```bash
git clone git@github.com:plantgenomicslab/HMMerMe.git
```

If an SSH key is not yet set up, please generate one. Navigate to the 'HMMerMe' folder:

```bash
cd HMMerMe
```

## Executing 'HMMerMe'

The main script is `run.py`, which requires two arguments: `--input` for the input directory and `--db` for the database directory. Optional arguments include `--CPU` to specify the number of cores (default is 2) and `--visualization` to generate data visualizations.

Execute the script as follows:

```bash
python run.py --input Input/ --db Database/ [OPTIONS --CPU, --visualization]
```

For data visualization, include the `--visualization` flag:

```bash
python run.py --input Input/ --db Database/ --visualization
```

### Input and Output Management

Expected input and output formats are detailed below, demonstrating the workflow from input fasta files and HMM profiles to the generation of various output files, including homology results, domain-specific lists, and visualizations of aligned sequences.

## Overview

This section elucidates how `.hmm` profile files in the `Database` directory and species fasta files in the `Input` directory work together within `run.py` to identify homologous sequences. It includes a detailed explanation of the output files generated post-analysis, their formats, and the significance of each file type in the context of HMM analysis.

### Custom Database Creation

For custom database creation, follow the steps below using `muscle` for alignment, `trimal` for trimming alignments based on conservation and gap thresholds, `esl-reformat` for format conversion, and `hmmbuild` to construct the HMM profile:

```bash
muscle -in {Homology_domain}.fasta -out {Homology_domain}.aln -clw
trimal -in {Homology_domain}.aln -out {Homology_domain}_trimmed.aln -gt 0.50 -cons 60
esl-reformat stockholm {Homology_domain}_trimmed.aln > {Homology_domain}.sto
hmmbuild {Homology_domain}.hmm {Homology_domain}.sto
```

These instructions guide you through setting up your analysis environment, executing the HMMerMe tool, and understanding the input and output files crucial for HMM analysis, offering a comprehensive toolkit for homology detection.
