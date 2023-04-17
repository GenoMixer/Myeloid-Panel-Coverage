# Myeloid-Panel-Coverage
Workflow for the Illumina TruSight Myeloid target-based coverage analysis using Snakemake.


## Requirement
- Linux (tested on Debian v9.13)
- snakemake v7.25.0   
- python v3.11.3 
- mosdepth v0.3.3
- bedtools v2.28.0
- r-xlsx v0.6.5


## General Information
For each Illumina TruSight Myeloid amplicon-based sequencing run, a target-based coverage analysis is performed to identify the mean coverage, regions and bases failed to reach the target region coverage. For this, a coverage threshold of 500x and exon/intron boundaries of +/-20 bp based on the coding sequences of the target regions were defined. This snakemake workflow will look for the samplesheet of the sequencing run and the "Alignment_1" directory with the output data from the Illumina Pisces pipeline and will serve as input for the coverage analysis. The pipeline should be executed from the corresponding flowcell directory and the output data will be stored in a toplevel directory called "coverage". The main results of the coverage analyses are stored in a sample-wise manner in excel-based worksheets. 


## Usage
1. Change directory to the MiSeq flowcell folder of the current sequencing run for which the coverage analysis should be run:

```bash
cd "/mnt/nas-5268189/ifh-rechenzentrum1/illumina/MiSeqOutput/230413"
```

2. Create a samplesheet with the corresponding samples and targeted genes as depicted here. First column with sample ids should be tab seperated and the corresponding genes comma seperated. The file should be named as following "Run_" + "Date" + "_Genliste.txt" and stored in the toplevel directory of the flowcell folder:

```bash
57204	ASXL1,TP53,RUNX1,IDH1,IDH2
57221	ASXL1,EZH2,TP53
57242	ASXL1,EZH2,TP53
57245	ASXL1,TP53,RUNX1,IDH1,IDH2
57268	ASXL1,EZH2,TP53,RUNX1,IDH1,IDH2
57280	ASXL1,EZH2,TP53,RUNX1,IDH1,IDH2
57296	ASXL1,EZH2,TP53
57307	ASXL1,TP53,RUNX1,IDH1,IDH2
```

3. Clone the repository into the runfolder:

    ```bash
    git clone "https://github.com/GenoMixer/Myeloid-Panel-Coverage.git"
    ```

4. Activate the snakemake environment with the relevant tools, i.e. mosdepth, bedtools and r-xlsx:

    ```bash
    conda activate "/mnt/nas-5268189/ifh-rechenzentrum1/bioinformatik/conda/envs/coverage"
    ```

5. Start a dry run (-n) and estimate the numbers of jobs(-j) provided by snakemake:

    ```bash
    snakemake -j 48 -n
    ```
