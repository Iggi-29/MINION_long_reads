##############
### FLAMES ###
##############

### align 
rule minimap_for_flames:
    input:
        fastq_filtered = base_path + "/flitered_fastq_files/{sample}/{sample}.fastq.gz"
    output:
        bam_flames_sorted = base_path + "/results/flames/{sample}_sorted.bam",
        bai_flames_sorted = base_path + "/results/flames/{sample}_sorted.bam.bai"
    params:
        genome = genome
    threads:
        10
    log:
        terminal_log = snake_files + "/log/{sample}_flames_alignment.log", 
        snakemake_log = snake_files + "/snakemake_log/{sample}_flames_alignment.log"
    benchmark:
        snake_files + "/benchmark/{sample}_flames_alignment.bmk"
    shell:"""
    minimap2 -t {threads} -ax splice -k14 --secondary=no --splice-flank=no {params.genome} {input.fastq_filtered} | \
    samtools view -bS -@ 4 -m 2G - | \
    samtools sort -@ {threads} -o {output.bam_flames_sorted}
    samtools index {output.bam_flames_sorted} {output.bai_flames_sorted}
    """

### flames
rule flames:
    input:
        bam_flames_sorted = base_path + "/results/flames/{sample}_sorted.bam",
        bai_flames_sorted = base_path + "/results/flames/{sample}_sorted.bam.bai",
        gff3_filtered = base_path + "/ref_files/filtered_MANE.gff3"
    output:
        FLAMES_done = base_path + "/results/flames/{sample}/FLAMES_done.txt"
    params:
        flames_script = flames_script, ### OJO
        genome        = genome,
        out_dir       = base_path + "/results/flames/{sample}",
        flames_config = flames_config, ### OJO
        minimap2_dir  = "/imppc/labs/eclab/ijarne/miniconda3/envs/FLAMES/bin/", ### OJO - NO CAL CANVIAR RES
        input_fasta = base_path + "/flitered_fastq_files/{sample}/{sample}.fastq.gz", 
        moved_input_fasta = base_path + "/results/flames/{sample}/merged.fastq.gz"
    threads: 
        10
    conda:
        flames_env
    log:
        terminal_log = snake_files + "/log/{sample}_flames.log", 
        snakemake_log = snake_files + "/snakemake_log/{sample}_flames.log"
    benchmark:
        snake_files + "/Run/benchmark/{sample}_flames.bmk"
    shell:"""
    ln -sf {params.input_fasta} {params.moved_input_fasta}
    
    cd {params.flames_script}

    ./python/bulk_long_pipeline.py \
    -a {input.gff3_filtered} \
    --genomefa {params.genome} \
    --outdir {params.out_dir} \
    --config_file {params.flames_config} \
    --minimap2_dir {params.minimap2_dir} \
    --inbam {input.bam_flames_sorted}

    echo "FLAMES has been performed on sample {wildcards.sample}" > {output.FLAMES_done}
    """

rule trigger_flames:
    input:
        expand(base_path + "/results/flames/{sample}/FLAMES_done.txt",
        base_path = base_path, sample = samples),
