#############
### FLAIR ###
#############

### align 
rule FLAIR_alignment:
    input:
        fastq_filtered = base_path + "/flitered_fastq_files/{sample}/{sample}.fastq.gz"
    output:
        # flair_aligned = base_path + "/results/alignments/flair/{sample}"
        done_step = base_path + "/results/alignments/flair/{sample}_flair.aligned_done.txt"
    params:
        flair_final =  base_path + "/results/alignments/flair/{sample}.flair.aligned",
        genome = genome,
    threads:
        9
    conda:
        flair_env
    log:
        terminal_log = snake_files + "/log/{sample}_flair_alignment.log", 
        snakemake_log = snake_files + "/snakemake_log/{sample}_flair_alignment.log"
    benchmark: 
        snake_files + "/benchmark/{sample}_flair_alignment.bmk"
    shell:"""
        flair align -g {params.genome} \
        -r {input.fastq_filtered} --output {params.flair_final} \
        --threads {threads}

        touch {output.done_step}
        """

### correction
rule FLAIR_correction:
    input:
        done_step = base_path + "/results/alignments/flair/{sample}_flair.aligned_done.txt"
    output:
        # flair_aligned = base_path + "/results/alignments/flair/{sample}"
        done_step = base_path + "/results/flair/correction/{sample}_flair.correction_done.txt"
    params:
        flair_final =  base_path + "/results/flair/correction/{sample}.flair.corrected",
        flair_aligned_bed = base_path + "/results/alignments/flair/{sample}.flair.aligned.bed",
        genome = genome,
        gtf_annotation = gtf_annotation
    threads:
        9
    conda:
        flair_env
    log:
        terminal_log = snake_files + "/log/{sample}_flair_correction.log", 
        snakemake_log = snake_files + "/snakemake_log/{sample}_flair_correction.log"
    benchmark: 
        snake_files + "/benchmark/{sample}_flair_correction.bmk"
    shell:"""
        flair correct -g {params.genome} \
        -f {params.gtf_annotation} \
        -q {params.flair_aligned_bed} \
        --output {params.flair_final} \
        --threads {threads}

        touch {output.done_step}
        """

### collapsing
rule FLAIR_collapse:
    input:
        done_step = base_path + "/results/flair/correction/{sample}_flair.correction_done.txt"
    output:
        # flair_aligned = base_path + "/results/alignments/flair/{sample}"
        done_step = base_path + "/results/flair/collapsed/{sample}_flair.collapsed_done.txt"
    params:
        flair_final = base_path + "/results/flair/collapsed/{sample}.flair.collapsed.bed",
        flair_corrected_bed = base_path + "/results/flair/correction/{sample}.flair.corrected_all_corrected.bed",
        fastq_filtered = base_path + "/flitered_fastq_files/{sample}/{sample}.fastq.gz",
        genome = genome,
        gtf_annotation = gtf_annotation
    threads:
        9
    conda:
        flair_env
    log:
        terminal_log = snake_files + "/log/{sample}_flair_collapsing.log", 
        snakemake_log = snake_files + "/snakemake_log/{sample}_flair_collapsing.log"
    benchmark: 
        snake_files + "/benchmark/{sample}_flair_collapsing.bmk"
    shell:"""
        flair collapse -g {params.genome} \
        -f {params.gtf_annotation} \
        -q {params.flair_corrected_bed} \
        -r {params.fastq_filtered} \
        --output {params.flair_final} \
        --threads {threads}

        touch {output.done_step}
        """

### quantify
rule FLAIR_quantify:
    input:
        done_step = expand(base_path + "/results/flair/collapsed/{sample}_flair.collapsed_done.txt",
        base_path = base_path, sample = samples)
    output:
        # flair_aligned = base_path + "/results/alignments/flair/{sample}"
        done_step = base_path + "/results/flair/quantify/{sample}_flair.quantify_done.txt"
    params:
        flair_final = base_path + "/results/flair/quantify/{sample}.flair.collapsed.bed",
        isoforms_bed_file = base_path + "/results/flair/collapsed/{sample}.flair.collapsed.bed.isoforms.fa",
        manifest_file =  manifest_file_lrs4,
    threads:
        20
    conda:
        flair_env
    log:
        terminal_log = snake_files + "/log/{sample}_flair_quantify.log", 
        snakemake_log = snake_files + "/snakemake_log/{sample}_flair_quantify.log"
    benchmark: 
        snake_files + "/benchmark/{sample}_flair_quantify.bmk"
    shell:"""
        flair quantify --sample_id_only -r {params.manifest_file} \
        -i {params.isoforms_bed_file} \
        --output {params.flair_final} \
        --threads {threads}

        touch {output.done_step}
        """

rule FLAIR_quantify2:
    input:
        done_step = expand(base_path + "/results/flair/collapsed/{sample}_flair.collapsed_done.txt",
        base_path = base_path, sample = samples)
    output:
        # flair_aligned = base_path + "/results/alignments/flair/{sample}"
        done_step = base_path + "/results/flair/quantify/{sample}_flair.quantify2_done.txt"
    params:
        flair_final = base_path + "/results/flair/quantify/{sample}.flair.collapsed.bed",
        isoforms_bed_file = base_path + "/results/flair/collapsed/{sample}.flair.collapsed.bed.isoforms.fa",
        manifest_file5 =  manifest_file_lrs5,
        manifest_file6 =  manifest_file_lrs6,
        manifest_file7 =  manifest_file_lrs7,
        manifest_file8 =  manifest_file_lrs8,
    threads:
        20
    conda:
        flair_env
    log:
        terminal_log = snake_files + "/log/{sample}_flair_quantify.log", 
        snakemake_log = snake_files + "/snakemake_log/{sample}_flair_quantify.log"
    benchmark: 
        snake_files + "/benchmark/{sample}_flair_quantify.bmk"
    shell:"""
        flair quantify --sample_id_only -r {params.manifest_file5} \
        -i {params.isoforms_bed_file} \
        --output {params.flair_final} \
        --threads {threads}

        flair quantify --sample_id_only -r {params.manifest_file6} \
        -i {params.isoforms_bed_file} \
        --output {params.flair_final} \
        --threads {threads}

        flair quantify --sample_id_only -r {params.manifest_file7} \
        -i {params.isoforms_bed_file} \
        --output {params.flair_final} \
        --threads {threads}

        flair quantify --sample_id_only -r {params.manifest_file8} \
        -i {params.isoforms_bed_file} \
        --output {params.flair_final} \
        --threads {threads}

        touch {output.done_step}
        """


rule trigger_flair:
    input:
        expand(base_path + "/results/alignments/flair/{sample}_flair.aligned_done.txt",
        base_path = base_path, sample = samples),
        expand(base_path + "/results/flair/correction/{sample}_flair.correction_done.txt",
        base_path = base_path, sample = samples),
        expand(base_path + "/results/flair/collapsed/{sample}_flair.collapsed_done.txt",
        base_path = base_path, sample = samples),
        expand(base_path + "/results/flair/quantify/{sample}_flair.quantify_done.txt",
        base_path = base_path, sample = samples_lrs4),
        expand(base_path + "/results/flair/quantify/{sample}_flair.quantify2_done.txt",
        base_path = base_path, sample = samples_lrs5 + samples_lrs6 + samples_lrs7 + samples_lrs8)


