##################################################
#### Snake file to Run Nanopore data analysis #### 
##################################################

import pandas as pd
from shutil import copyfile
import os, glob, re,sys

configfile:
    "config/config.yaml"

samples = [s for s in config["samples"]]
base_path = "/imppc/labs/eclab/ijarne/0_Recerca/4_MINION_better"
concatenated_path = "/imppc/labs/eclab/ijarne/0_Recerca/4_MINION_better/fastq_concat"
genome = "/imppc/labs/eclab/Resources/MINION_ddbb/flair_execution/GRCh38.primary_assembly.genome.fa"
transcriptome = "/imppc/labs/eclab/Resources/MINION_ddbb/flair_execution/gencode.v48.pc_transcripts.fa"
gtf_file = "/imppc/labs/eclab/Resources/MINION_ddbb/flair_execution/gencode.v46.basic.annotation.gtf"
genes_to_look_at = ["LZTR1","NF1","NF2","PKD2","PKD1","SMARCE1","SMARCB1","TSC1","TSC2","SPRED1"]
# FAMES_sample_trial = ["barcode03"]
FLAMES_sample_trial = ["barcode15","barcode16"]

rule all:
    input:
        ####### Fastp
        expand(base_path + "/results/1_fastp/{sample}/{sample}.fastq.gz",
        base_path = base_path, sample = samples),
        ####### BAMSLAM
        expand(base_path + "/results/3_Alignment/transcriptome/{sample}.bam",
        base_path = base_path, sample = samples),
        expand(base_path + "/results/1_bam_slam/{sample}/bam_slam_done.txt",
        base_path = base_path, sample = samples),
        ####### BAMSLAM - per gene
        expand(base_path + "/results/1_bam_slam_per_gene/{sample}/{gene}/bam_slam_done.txt",
        base_path = base_path, sample = samples, gene = genes_to_look_at),
        ####### Coverage report
        # base_path + "/results/1_coverage_report/cov_report.tsv",
        base_path + "/results/1_coverage_report/Heatmaps_done.txt",
        base_path + "/results/1_coverage_report/Barplots_done.txt",
        ####### GFF3 modification - prepare its execution
        # expand("/imppc/labs/eclab/Resources/MINION_ddbb/flair_execution/per_gene/{gene}.gff3",
        # gene = genes_to_look_at),
        ####### FLAMES
        # expand(base_path + "/results/3_Alignment/genome_FLAMES/{sample}_sorted.bam",
        # base_path = base_path, sample = FLAMES_sample_trial),
        # expand(base_path + "/results/4_FLAMES/{sample}/FLAMES_done.txt",
        # base_path = base_path, sample = FAMES_sample_trial),
        expand(base_path + "/results/4_FLAMES/{sample}/FLAMES_done.txt",
        base_path = base_path, sample = FLAMES_sample_trial),
        ####### BAMBU
        expand(base_path + "/results/3_Alignment/genome/{sample}_sorted.bam",
        base_path = base_path, sample = samples),
        expand(base_path + "/results/3_Alignment/genome/{sample}_sorted.bam.bai",
        base_path = base_path, sample = samples),


######################################################################################
##                                   QC                                             ##
######################################################################################
########################
##      fastplong     ##
########################
rule fastp:
    input:
        concat_reads = lambda wc:config["samples"][wc.sample]["concat_reads"]
    output:
        json_out = base_path + "/results/1_fastp/{sample}/{sample}_report.json",
        html_out = base_path + "/results/1_fastp/{sample}/{sample}_report.html",
        fastq_filtered = base_path + "/results/1_fastp/{sample}/{sample}.fastq.gz"
    params:
        length_required = 500,
        q = 10,
        cut_front_mean_quality = 8,
        cut_tail_mean_quality = 8
    threads:
        6
    log:
        base_path + "/Run/logs/{sample}/fastp.log"
    benchmark:
        base_path + "/Run/benchmarks/{sample}/fastp.bmk"
    conda:
        "/imppc/labs/eclab/ijarne/0_Recerca/4_MINION_better/envs/fastp_qc.yaml"
    shell:"""
        fastplong --length_required {params.length_required} \
        --thread {threads} \
        --cut_front --cut_tail --cut_front_mean_quality {params.cut_front_mean_quality} \
        --cut_tail_mean_quality {params.cut_tail_mean_quality} \
        --in {input.concat_reads} --out {output.fastq_filtered} \
        --html {output.html_out} --json {output.json_out}
        """
##### just in case we have here code for the filtlong algorithm:
##### filtlong --min_length 500 \
##### --min_mean_q 10 \
##### /imppc/labs/eclab/ijarne/0_Recerca/4_MINION_better/concat_data/barcode03/barcode03.fastq.gz 
#####| gzip > /imppc/labs/eclab/ijarne/0_Recerca/4_MINION_better/results/1_filt_long/barcode03.fastq.gz

########################
## BamSlam whole data ##
########################
rule minimap_transcriptome:
    input:
        fastq_filtered = base_path + "/results/1_fastp/{sample}/{sample}.fastq.gz"
    output:
        bam_transcriptome = base_path + "/results/3_Alignment/transcriptome/{sample}.bam"
    params:
        transcriptome = transcriptome
    threads:
        6
    log:
        base_path + "/Run/logs/{sample}/minimap_transcriptome.log"
    benchmark:
        base_path + "/Run/benchmarks/{sample}/minimap_transcriptome.bmk"
    shell:"""
    minimap2 -ax map-ont --sam-hit-only {params.transcriptome} {input.fastq_filtered} | samtools view -bh >  {output.bam_transcriptome}
    """

rule bam_slam:
    input:
        bam_transcriptome = base_path + "/results/3_Alignment/transcriptome/{sample}.bam"
    output:
        bam_slam_done = base_path + "/results/1_bam_slam/{sample}/bam_slam_done.txt"
    params:
        bam_slam_script = base_path + "/scripts/BamSlam.R",
        bam_slam_out_name = base_path + "/results/1_bam_slam/{sample}/"
    threads:
        1
    conda:
        "/imppc/labs/eclab/ijarne/0_Recerca/4_MINION_better/envs/g_alignments.yaml"
    log:
        base_path + "/Run/logs/{sample}/bam_slam.log"
    benchmark:
        base_path + "/Run/benchmarks/{sample}/bam_slam.bmk"
    shell:"""
    Rscript {params.bam_slam_script} cdna {input.bam_transcriptome} {params.bam_slam_out_name} {output.bam_slam_done}
    """
#### in case: 
#### Rscript 
# ./scripts/BamSlam.R cdna /imppc/labs/eclab/ijarne/0_Recerca/4_MINION_better/results/3_Alignment/transcrptome/barcode05.bam
# /imppc/labs/eclab/ijarne/0_Recerca/4_MINION_better/results/3_Alignment/transcrptome/barcode05

###########################
## BamSlam for each gene ##
###########################
rule gene_minimap_transcriptome:
    input:
        fastq_filtered = base_path + "/results/1_fastp/{sample}/{sample}.fastq.gz"
    output:
        bam_transcriptome = temp(base_path + "/results/3_Alignment/transcriptome/{sample}_{gene}.bam")
    params:
        transcriptome = "/imppc/labs/eclab/Resources/MINION_ddbb/flair_execution/per_gene/{gene}.fasta"
    threads:
        2
    log:
        base_path + "/Run/logs/{sample}/minimap_transcriptome_{gene}.log"
    benchmark:
        base_path + "/Run/benchmarks/{sample}/minimap_transcriptome_{gene}.bmk"
    shell:"""
    minimap2 -ax map-ont --sam-hit-only {params.transcriptome} {input.fastq_filtered} | samtools view -bh >  {output.bam_transcriptome}
    """

rule gene_bam_slam:
    input:
        bam_transcriptome = base_path + "/results/3_Alignment/transcriptome/{sample}_{gene}.bam"
    output:
        bam_slam_done = base_path + "/results/1_bam_slam_per_gene/{sample}/{gene}/bam_slam_done.txt"
    params:
        bam_slam_script = base_path + "/scripts/BamSlam.R",
        bam_slam_out_name = base_path + "/results/1_bam_slam_per_gene/{sample}/{gene}/"
    threads:
        1
    conda:
        "/imppc/labs/eclab/ijarne/0_Recerca/4_MINION_better/envs/g_alignments.yaml"
    log:
        base_path + "/Run/logs/{sample}/{gene}_bam_slam.log"
    benchmark:
        base_path + "/Run/benchmarks/{sample}/{gene}_bam_slam.bmk"
    shell:"""
    Rscript {params.bam_slam_script} cdna {input.bam_transcriptome} {params.bam_slam_out_name} {output.bam_slam_done}
    """

####################
##    Coverage    ##
####################
# rule general_align:
#     input:
#         fastq_filtered = base_path + "/results/1_fastp/{sample}/{sample}.fastq.gz"
#     output:
#         bam_genome_sorted= temp(base_path + "/results/3_Alignment/genome_coverage/{sample}_sorted.bam"),
#         bai_genome_sorted = temp(base_path + "/results/3_Alignment/genome_coverage/{sample}_sorted.bam.bai"),
#     params:
#         genome = "/imppc/labs/eclab/Resources/MINION_ddbb/flair_execution/GRCh38.primary_assembly.genome.fa"
#     threads:
#         6
#     log:
#         base_path + "/Run/logs/{sample}/minimap_genome_coverage.log"
#     benchmark:
#         base_path + "/Run/benchmarks/{sample}/minimap_genome_coverage.bmk"
#     shell:"""
#     minimap2 -ax splice --secondary=no --sam-hit-only {params.genome} {input.fastq_filtered} | \
#     samtools view -bh - | \
#     samtools sort -@ {threads} -o {output.bam_genome_sorted}
#     samtools index {output.bam_genome_sorted} {output.bai_genome_sorted}
#     """
rule exon_coverage_report:
    input:
        bam_genome_sorted = expand(base_path + "/results/3_Alignment/genome/{sample}_sorted.bam",
        base_path = base_path, sample = samples),
        bai_genome_sorted = expand(base_path + "/results/3_Alignment/genome/{sample}_sorted.bam.bai",
        base_path = base_path, sample = samples),
    output:
        final_coverage_report = base_path + "/results/1_coverage_report/cov_report.tsv"
    params:
        bed_files_folder = base_path + "/bed_files",
        sam_files_folder = base_path + "/results/3_Alignment/genome",
        coverage_folder = base_path + "/results/1_coverage",
        scripts_place = base_path + "/scripts"
    threads:
        2
    log:
        base_path + "/Run/logs/Run/exon_coverage_report.log"
    benchmark:
        base_path + "/Run/benchmarks/Run/exon_coverage_report.bmk"
    conda:
        "/imppc/labs/eclab/ijarne/0_Recerca/4_MINION_better/envs/bed_generation_R.yaml"
    script:"{params.scripts_place}/1_exon_coverage.R"


rule exon_coverage_heatmaps:
    input:
        final_coverage_report = base_path + "/results/1_coverage_report/cov_report.tsv"
    output:
        heatmaps_done = base_path + "/results/1_coverage_report/Heatmaps_done.txt"
    params:
        scripts_place = base_path + "/scripts",
        coverage_folder = base_path + "/results/1_coverage_report",
    threads:
        2
    log:
        base_path + "/Run/logs/Run/exon_coverage_HM.log"
    benchmark:
        base_path + "/Run/benchmarks/Run/exon_coverage_HM.bmk"
    conda:
        "/imppc/labs/eclab/ijarne/0_Recerca/4_MINION_better/envs/g_alignments.yaml"
    script:"{params.scripts_place}/2_exon_coverage_plots.R"

rule exon_coverage_barplots:
    input:
        final_coverage_report = base_path + "/results/1_coverage_report/cov_report.tsv"
    output:
        bars_done = base_path + "/results/1_coverage_report/Barplots_done.txt"
    params:
        scripts_place = base_path + "/scripts",
        coverage_folder = base_path + "/results/1_coverage_report",
    threads:
        2
    log:
        base_path + "/Run/logs/Run/exon_coverage_barplots.log"
    benchmark:
        base_path + "/Run/benchmarks/Run/exon_coverage_barplots.bmk"
    conda:
        "/imppc/labs/eclab/ijarne/0_Recerca/4_MINION_better/envs/g_alignments.yaml"
    script:"{params.scripts_place}/3_exon_coverage_plots.R"


###################################################
## BAMBU pipeline, alignment and bambu execution ##
###################################################
rule minimap2_for_bambu:
    input:
        fastq_filtered = base_path + "/results/1_fastp/{sample}/{sample}.fastq.gz"
    output:
        # bam_genome = base_path + "/results/3_Alignment/genome/{sample}.bam",
        bam_genome_sorted = base_path + "/results/3_Alignment/genome/{sample}_sorted.bam",
        bai_genome_sorted = base_path + "/results/3_Alignment/genome/{sample}_sorted.bam.bai"
    params:
        genome = genome
    threads:
        6
    log:
        base_path + "/Run/logs/{sample}/minimap_genome.log"
    benchmark:
        base_path + "/Run/benchmarks/{sample}/minimap_genome.bmk"
    shell:"""
    minimap2 -t {threads} -ax splice -uf --MD --secondary=no {params.genome} {input.fastq_filtered} | \
    samtools view -bh - | \
    samtools sort -@ {threads} -o {output.bam_genome_sorted}
    samtools index {output.bam_genome_sorted} {output.bai_genome_sorted}
    """

# minimap2 -ax splice -u --MD --secondary=no \
#   /imppc/labs/eclab/Resources/MINION_ddbb/flair_execution/GRCh38.primary_assembly.genome.fa \
#   /imppc/labs/eclab/ijarne/0_Recerca/4_MINION_better/results/1_fastp/barcode03/barcode03.fastq.gz | \
#   samtools view -bh - | \
#   samtools sort -o /imppc/labs/eclab/ijarne/0_Recerca/4_MINION_better/results/3_Alignment/genome/barcode03_sorted.bam

# rule bambu_isforms:
#     input: 
#         ""
#     output: 
#         ""
#     params: 
#         ""
#     threads:
#         6
#     log:
#         base_path + "/Run/logs/{sample}/bambu_isforms.log"
#     benchmark:
#         base_path + "/Run/benchmarks/{sample}/bambu_isforms.bmk"
#     shell:
#         ""
###########################
## GFF3 filter with AGAT ##
###########################
rule GFF3_filter:
    input:
        txt_of_filter = "/imppc/labs/eclab/ijarne/0_Recerca/4_MINION_better/config/keep_this.txt",
        gff3 = "/imppc/labs/eclab/Resources/MINION_ddbb/flair_execution/gencode.v46.annotation.gff3",
    output:
        gff3_filtered = "/imppc/labs/eclab/Resources/MINION_ddbb/flair_execution/filtered.gff3"
    conda:
        "/imppc/labs/eclab/ijarne/0_Recerca/4_MINION_better/envs/AGAT.yaml"
    log:
        base_path + "/Run/logs/Run/gff_filter.log"
    benchmark:
        base_path + "/Run/benchmarks/Run/gff_filter.bmk"
    threads:
        1
    shell:"""
    agat_sp_filter_feature_from_keep_list.pl -f {input.gff3} -a gene_name -kl {input.txt_of_filter} -o {output.gff3_filtered}
    """
        
############
## FLAMES ##
############
rule minimap2_for_flames:
    input:
        fastq_filtered = base_path + "/results/1_fastp/{sample}/{sample}.fastq.gz"
    output:
        bam_genome_sorted = base_path + "/results/3_Alignment/genome_FLAMES/{sample}_sorted.bam",
        bai_genome_sorted = base_path + "/results/3_Alignment/genome_FLAMES/{sample}_sorted.bam.bai"
    params:
        genome = genome
    threads:
        10
    log:
        base_path + "/Run/logs/{sample}/flames_alignment.log"
    benchmark:
        base_path + "/Run/benchmarks/{sample}/flames_alignment.bmk"
    shell:"""
    minimap2 -t {threads} -ax splice -k14 --secondary=no --splice-flank=no {params.genome} {input.fastq_filtered} | \
    samtools view -bS -@ 4 -m 2G - | \
    samtools sort -@ {threads} -o {output.bam_genome_sorted}
    samtools index {output.bam_genome_sorted} {output.bai_genome_sorted}
    """

rule flames:
    input:
        bam_genome_sorted = base_path + "/results/3_Alignment/genome_FLAMES/{sample}_sorted.bam",
        gff3_filtered = "/imppc/labs/eclab/Resources/MINION_ddbb/flair_execution/filtered.gff3"
    output:
        FLAMES_done = base_path + "/results/4_FLAMES/{sample}/FLAMES_done.txt"
    params:
        FLAMES_script = "/imppc/labs/eclab/ijarne/0_Recerca/4_MINION_better/scripts/FLAMES",
        genome        = genome,
        out_dir       = base_path + "/results/4_FLAMES/{sample}",
        FLAMES_config = "/imppc/labs/eclab/ijarne/0_Recerca/4_MINION_better/config/FLAMES_config.json",
        minimap2_dir  = "/imppc/labs/eclab/ijarne/miniconda3/envs/FLAMES/bin/",
        input_fasta = base_path + "/results/1_fastp/{sample}/{sample}.fastq.gz",
        moved_input_fasta = base_path + "/results/4_FLAMES/{sample}/merged.fastq.gz"
    threads: 10
    conda:
        "/imppc/labs/eclab/ijarne/0_Recerca/4_MINION_better/envs/FLAMES.yaml"
    log:
        base_path + "/Run/logs/{sample}/flames.log"
    benchmark:
        base_path + "/Run/benchmarks/{sample}/flames.bmk"
    shell:"""
    ln -sf {params.input_fasta} {params.moved_input_fasta}
    
    cd {params.FLAMES_script}

    ./python/bulk_long_pipeline.py \
    -a {input.gff3_filtered} \
    --genomefa {params.genome} \
    --outdir {params.out_dir} \
    --config_file {params.FLAMES_config} \
    --minimap2_dir {params.minimap2_dir} \
    --inbam {input.bam_genome_sorted}

    echo "FLAMES has been performed on sample {wildcards.sample}" > {output.FLAMES_done}
    """
    
# rule flames:
#     input:
#         fastq_filtered = base_path + "/results/1_fastp/{sample}/{sample}.fastq.gz"
#     output:
#         FLAMES_done = base_path + "/results/4_FLAMES/{sample}/FLAMES_done.txt"
#     params:
#         FLAMES_script = "/imppc/labs/eclab/ijarne/0_Recerca/4_MINION_better/scripts/FLAMES",
#         gff3 = "/imppc/labs/eclab/Resources/MINION_ddbb/flair_execution/gencode.v46.annotation.gff3",
#         genome = "/imppc/labs/eclab/Resources/MINION_ddbb/flair_execution/GRCh38.primary_assembly.genome.fa",
#         out_dir = base_path + "/results/4_FLAMES/{sample}",
#         FLAMES_config = "/imppc/labs/eclab/ijarne/0_Recerca/4_MINION_better/config/FLAMES_config.json",
#         minimap2_dir = "/imppc/labs/eclab/ijarne/miniconda3/envs/FLAMES/bin/",
#         fastq_dir = base_path + "/results/1_fastp/{sample}"
#     threads:
#         6
#     conda:
#         "/imppc/labs/eclab/ijarne/0_Recerca/4_MINION_better/envs/FLAMES.yaml"
#     log:
#         base_path + "/Run/logs/{sample}/flames.log"
#     benchmark:
#         base_path + "/Run/benchmarks/{sample}/flames.bmk"
#     shell:"""
#         cd {params.FLAMES_script}
#         ./python/bulk_long_pipeline.py -a {params.gff3} \
#         --genomefa {params.genome} \
#         --outdir {params.out_dir} \
#         --config_file {params.FLAMES_config} \
#         --minimap2_dir  {params.minimap2_dir}\
#         --fq_dir {params.fastq_dir}
#         echo "FLAMES has been performed on this sample {sample}" > {output.FLAMES_done}
#         """

