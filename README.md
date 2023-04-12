# Myeloid-Panel-Coverage
Workflow for the Illumina TruSight Myeloid target-based coverage analysis using Snakemake.

## Requirement
- Linux (tested on Debian 9.13)
- snakemake v7.25.0   
- python 3.5+ 
- mosdepth v0.3.3
- bedtools v2.28.0
- r-xlsx v0.6.5


## General Information


When a NextSeq 550 sequencing run is completed the FASTQ files are located in "<run folder>\Alignment_1\<subfolder>\Fastq". For each sample the sequencer generates fastq files per lane and read orientation. The fastq files are stored initially in this format ("*_S[1-8]_L00[1-4]_R[1-2]_001.fastq.gz"). After merging the fastq files the Illumina-specific sample and lane information will be discarded. The output of the merged fastq files will be stored in a directory called "unaligned" with the corresponding QC and checksums. The resulting directory structure is highlighted here:
