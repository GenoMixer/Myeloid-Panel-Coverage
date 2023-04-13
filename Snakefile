# libraries
import os
import io
import glob
import numpy as np
import pandas as pd
import os.path as path

# config
#configfile: "config.yaml"
WORK  = os.getcwd()
head, tail = os.path.split(WORK)
RUN = tail.split("_")[0]
print(RUN)

# skip until match 
class SkipUntilMatchWrapper(io.TextIOWrapper):
    def __init__(self, f, matcher, include_matching=False):
        super().__init__(f, line_buffering=True)
        self.f = f
        self.matcher = matcher
        self.include_matching = include_matching
        self.has_matched = False

    def read(self, size=None):
        while not self.has_matched:
            line = self.readline()
            if self.matcher(line):
                self.has_matched = True
                if self.include_matching:
                    return line
        return super().read(size)

# read samples
with open('SampleSheet.csv', 'rb') as f_orig:
    with SkipUntilMatchWrapper(f_orig, lambda s: '[Data]' in s, include_matching=True) as f:
        samples=pd.read_csv(f, skiprows=1).set_index("Sample_ID")
print(samples)

# get bam
def get_bam(wildcards):
        return sorted(glob.glob("Alignment_1/*/" + wildcards.sample + "_S[1-8]" + ".bam"))

rule all:
    input:
         expand("coverage/{sample}.mosdepth.global.dist.txt", sample=samples.index),
         expand("coverage/{sample}.per-base.bed.gz", sample=samples.index),
         expand("coverage/{sample}.per-base.bed", sample=samples.index),
         expand("coverage/{sample}.mosdepth.summary.txt", sample=samples.index),
         expand("coverage/{sample}.thresholds.bed.gz", sample=samples.index),
         expand("coverage/{sample}.thresholds.bed", sample=samples.index),
         expand("coverage/{sample}.regions.bed.gz", sample=samples.index),
         expand("coverage/{sample}.regions.bed", sample=samples.index),
         expand("coverage/{sample}.targets_pct_gt_500.bed", sample=samples.index),
         expand("coverage/{sample}.targets_not_100pct_gt_500.bed", sample=samples.index),
         expand("coverage/{sample}.targets_not_100pct_gt_500.txt", sample=samples.index),
         expand("coverage/{sample}.bases_lt_500_1based_annotated_regions.txt", sample=samples.index),
         expand("coverage/{sample}.genelist", sample=samples.index),
         expand("coverage/{sample}.genelist_cds_lengths.txt", sample=samples.index),
         expand("coverage/{sample}.regions_meancov_genelist.bed", sample=samples.index),
         expand("coverage/{sample}.bases_lt_500_1based_annotated_genelist.txt", sample=samples.index),
         expand("coverage/{sample}.bases_lt_500_1based_annotated_regions_genelist.txt", sample=samples.index),
         expand("coverage/{sample}.targets_not_100pct_gt_500_genelist.txt", sample=samples.index),
         expand("coverage/{sample}.coverage_data_compilation.xlsx", sample=samples.index)


rule mosdepth:
    input:
        bam=get_bam,
        bed="/mnt/ngs-resources/transcript_annotations/refseq_hg19/TruSight-Myeloid-Amplicon-Panel/TruSight-Myeloid-Amplicon-Panel.hg19.selected_transcripts_cds.sorted.int20.intersect_amplicons.bed"
    output:
        "coverage/{sample}.mosdepth.global.dist.txt",
        "coverage/{sample}.per-base.bed.gz",
        "coverage/{sample}.mosdepth.summary.txt",
        "coverage/{sample}.thresholds.bed.gz",
        "coverage/{sample}.regions.bed.gz"
    params:
        extra="--fast-mode --mapq 20 --thresholds 0,1,500,1000",
        prefix="{sample}"
    threads: 4
    shell: 
        """
        mosdepth --by {input.bed} --threads {threads} {params.extra} coverage/{params.prefix} {input.bam}
        """

rule pigz1:
    input: "coverage/{sample}.thresholds.bed.gz",
    output: "coverage/{sample}.thresholds.bed",
    shell: 
        """
        pigz -kdf {input} > {output}
        """


rule pigz2:
    input: "coverage/{sample}.per-base.bed.gz"
    output: "coverage/{sample}.per-base.bed"
    shell:
        """
        pigz -kdcf {input} > {output} # grep -v "HLA"
        """


rule pigz3:
    input: "coverage/{sample}.regions.bed.gz"
    output: "coverage/{sample}.regions.bed"
    shell:
        """
        pigz -kdf {input} > {output}
        """


rule create_bed:
    input: "coverage/{sample}.thresholds.bed"
    output: "coverage/{sample}.targets_pct_gt_500.bed"
    shell:
        """
        awk 'BEGIN {{FS="\t"; OFS="\t"}} NR>1 {{print $1,$2,$3,$4,($7/$5)*100}}' {input} > {output}
        """


rule create_bed_two:
    input: "coverage/{sample}.targets_pct_gt_500.bed"
    output: "coverage/{sample}.targets_not_100pct_gt_500.bed"
    shell:
        """
        awk 'BEGIN {{FS="\t"; OFS="\t"}} {{if ($5 < 100) print $0}}' {input} > {output}
        """


rule create_bed_r:
    input: "coverage/{sample}.targets_not_100pct_gt_500.bed"
    output: "coverage/{sample}.targets_not_100pct_gt_500.txt"
    shell: 
        """
        awk 'BEGIN {{FS="\t"; OFS="\t"}} {{print $1,$2+1,$3,$4,$5}}' {input} > {output}
        """


rule intersect:
    input:
        abed="coverage/{sample}.per-base.bed",
        bbed="coverage/{sample}.targets_not_100pct_gt_500.bed"
    output: "coverage/{sample}.bases_lt_500_1based.txt"
    shell: 
        """
        bedtools intersect \
           -a {input.abed} \
           -b {input.bbed} \
        | awk 'BEGIN {{FS="\t"; OFS="\t"}} {{if ($4 < 500) print $0}}' \
        | bedtools makewindows \
           -w 1 \
           -b - \
        | bedtools intersect \
           -a {input.abed} \
           -b - \
        | awk 'BEGIN {{FS="\t"; OFS="\t"}} {{print $1,$2+1,$3,$4}}' \
        | uniq > {output}
        """


rule intersect_two:
    input: 
        abed="coverage/{sample}.bases_lt_500_1based.txt",
        bbed="/mnt/ngs-resources/transcript_annotations/refseq_hg19/TruSight-Myeloid-Amplicon-Panel/TruSight-Myeloid-Amplicon-Panel.hg19.selected_transcripts_cds.sorted.int20.intersect_amplicons.bed.all_bases_annotated.txt"
    output: "coverage/{sample}.bases_lt_500_1based_annotated.txt"
    shell:
        """
        bedtools intersect \
           -a {input.abed} \
           -b {input.bbed} \
           -f 1 \
           -wb \
        | awk 'BEGIN {{
           FS="\t"
           OFS="\t"
           print "#chr","start","end","cov","gene","transcript","exon","coding_change","protein_consequence","region_type"
           }}
           {{
           sub("exonic","cds",$13)
           print $5,$6,$7,$4,$8,$9,$10,$11,$12,$13
           }}' - \
        | uniq > {output}
        """


rule sort:
    input: "coverage/{sample}.bases_lt_500_1based_annotated.txt"
    output: "coverage/{sample}.bases_lt_500_1based_annotated_regions.txt"
    shell:
        """
        sort -k1,1V -k2,2n -k3,3n \
           {input} \
        | bedtools merge \
           -i - \
           -c 5,6,7,8,8,9,9,10,10,4,4,4 -o distinct,distinct,distinct,first,last,first,last,first,last,mean,min,max \
        | awk 'BEGIN {{
           FS="\t"
           OFS="\t"
           print "#chr","first_base","last_base","gene","transcript","exon","first_coding_change","last_coding_change","first_protein_cons","last_protein_cons","region_length_bp","mean_cov","min_cov","max-cov"
        }}
        {{
          if ($3-$2 > 0) {{
             start_coord = $2+1
             end_coord = $3-1
        }}
          else {{
             start_coord = $2
             end_coord = $3
        }}
        exon = $6
        if ($exon ~ ".,") gsub(".,","",exon)
        first_prot = $9
        last_prot = $10
        if ($9 == "." && 11 != "cds") first_prot = $11
        if ($10 == "." && $12 != "cds") last_prot = $12
        print $1,start_coord,end_coord,$4,$5,exon,$7,$8,first_prot,last_prot,(end_coord-start_coord)+1,$13,$14,$15
        }}' - > {output}
        """


rule genelists:
    input: "Run_" + RUN + "_Genliste.txt"
    output: "coverage/{sample}.genelist"
    params: prefix="{sample}"
    shell:   
        """
        grep -w {params.prefix} {input} \
        | awk 'BEGIN {{FS="\t"; OFS="\t"}}
        {{
           n_genes=split($2,genes,",")
           for (i=1; i <= n_genes; i++ )
             print genes[i]
        }}' - > {output}
        """


rule cds_length:
    input:
        genes="coverage/{sample}.genelist",
        cds="/mnt/ngs-resources/transcript_annotations/refseq_hg19/TruSight-Myeloid-Amplicon-Panel/TruSight-Myeloid-Amplicon-Panel.hg19.genes_selected_transcripts_intersect_amplicons_cdslengths.bed"
    output: "coverage/{sample}.genelist_cds_lengths.txt"
    shell:
        """
        grep -F \
        -f {input.genes} \
        -w \
        {input.cds} \
        > {output}
        """


rule meancov:
    input:
        genes="coverage/{sample}.genelist",
        regions="coverage/{sample}.regions.bed"
    output: "coverage/{sample}.regions_meancov_genelist.bed"
    shell: 
        """
        grep -F \
        -f {input.genes} \
        {input.regions} \
        | awk 'BEGIN {{FS="\t"; OFS="\t"; print "Chr","Start","End","Target_Name","Mean_Coverage"}} {{print $1,$2+1,$3,$4,$5}}' - \
        > {output}
        """


rule basecov:
    input:
        genes="coverage/{sample}.genelist",
        regions="coverage/{sample}.bases_lt_500_1based_annotated.txt"
    output: "coverage/{sample}.bases_lt_500_1based_annotated_genelist.txt"
    shell:
        """
        set -e \
        grep -F \
        -f {input.genes} \
        {input.regions} \
        | awk 'BEGIN {{FS="\t"; OFS="\t"; print "#chr","start","end","cov","gene","transcript","exon","coding_change","protein_consequence","region_type"}} {{print $0}}' - \
        > {output}
        """


rule regions:
    input:
        genes="coverage/{sample}.genelist",
        regions="coverage/{sample}.bases_lt_500_1based_annotated_regions.txt"
    output: "coverage/{sample}.bases_lt_500_1based_annotated_regions_genelist.txt"
    shell:
        """
        set -e \
        grep -F \
        -f {input.genes} \
        {input.regions} \
        | awk 'BEGIN {{FS="\t"; OFS="\t"; print "#chr","first_base","last_base","gene","transcript","exon","first_coding_change","last_coding_change","first_protein_cons","last_protein_cons","region_length_bp","mean_cov","min_cov","max-cov"}} {{print $0}}' - \
        > {output}
        """


rule base:
    input:
        genes="coverage/{sample}.genelist",
        target="coverage/{sample}.targets_not_100pct_gt_500.txt"
    output: "coverage/{sample}.targets_not_100pct_gt_500_genelist.txt"
    shell: 
        """
        set -e \
        grep -F \
        -f {input.genes} \
        {input.target} \
        | awk 'BEGIN {{FS="\t"; OFS="\t"; print "Chr","Start","End","Target_Name","pct_bases_gt_30x_cov"}} {{print $0}}' - \
        > {output}
        """


rule rscript:
    input: 
        cds="coverage/{sample}.genelist_cds_lengths.txt",
        cov="coverage/{sample}.regions_meancov_genelist.bed",
        pct="coverage/{sample}.targets_not_100pct_gt_500_genelist.txt",
        regions="coverage/{sample}.bases_lt_500_1based_annotated_regions_genelist.txt",
        bases="coverage/{sample}.bases_lt_500_1based_annotated_genelist.txt"
    output: 
        out="coverage/{sample}.coverage_data_compilation.xlsx"
    params:  
       prefix="{sample}",
       run=RUN
#    conda: "xlsx.yml"
    script: 'R_coverage_data_compilation.R'

