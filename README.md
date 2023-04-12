# Myeloid-Panel-Coverage
Workflow for the Illumina TruSight Myeloid target-based coverage analysis using Snakemake.


## General Information
For each Illumina TruSight Myeloid Panel sequencing run on the MiSeq, a target-based coverage analysis is performed to identify the mean coverage, regions and bases failed to reach the target coverage. For this, a coverage threshold of 500x and exon/intron boundaries of +/-20bp based on the coding sequences of the target regions were defined. The snakemake workflow will be executed from the corresponding flowcell directory and the output data will be stored in a toplevel directory called "coverage". The main results of the coverage analyses are stored in sample-wise manner in excel-based worksheets. 


## Requirement
- Linux (tested on Debian v9.13)
- snakemake v7.25.0   
- python v3.5+ 
- mosdepth v0.3.3
- bedtools v2.28.0
- r-xlsx v0.6.5


## Usage

1. Clone the repository into the runfolder:

    ```bash
    git clone https://github.com/GenoMixer/Myeloid-Panel-Coverage.git
    ```

2. Activate the snakemake environment with the relevant tools, i.e. mosdepth, bedtools and r-xlsx:

    ```bash
    conda activate /mnt/nas-5268189/ifh-rechenzentrum1/bioinformatik/conda/envs/coverage
    ```

3. Start a dry run (-n) and estimate the numbers of jobs(-j) provided by snakemake:

    ```bash
    snakemake -j 48 -n
    ```
