import re
import subprocess
import os

#configfile: "misc/masterconfig2.yaml"

# Get raw_reads
#sample_ids, = glob_wildcards(config['raw_reads_path']+"/{sample}.R1.fastq.gz")
prefix = config['prefix']
maxthreads = snakemake.utils.available_cpu_count()

logs = config['base_log_outdir']

#print(sample_ids)

# one rule to rule them all
#rule all:
#    input:
#        expand(config['outdir']+"/{prefix}/fastp/{sample}.R1.fastq.gz",
#               sample=sample_ids, prefix=prefix),
#        expand(config['outdir']+"/{prefix}/summaries/fastp_summary.json", prefix=prefix)


rule run_fastp:
    input:
        r1 = config['raw_reads_path']+"/{sample}.R1.fastq.gz",
        r2 = config['raw_reads_path']+"/{sample}.R2.fastq.gz"
    output:
        r1_filt = config['outdir']+"/{prefix}/fastp/{sample}.R1.fastq.gz",
        r2_filt = config['outdir']+"/{prefix}/fastp/{sample}.R2.fastq.gz",
        json = config['outdir']+"/{prefix}/fastp/{sample}.fastp.json",
        html = config['outdir']+"/{prefix}/fastp/{sample}.fastp.html"
    log:
        config['base_log_outdir']+"/{prefix}/fastp/run/{sample}.log"
    threads:
        2
    conda:
        "config/fastp.yaml"
    shell:
        """
        fastp -i {input.r1} -I {input.r2} -o {output.r1_filt} -O {output.r2_filt} -j {output.json} -h {output.html}
        """

rule run_fastp_summary:
    input:
        expand(config['outdir']+"/{prefix}/fastp/{sample}.fastp.json", sample=sample_ids, prefix=prefix)
    output:
        config['outdir']+"/{prefix}/summaries/fastp_summary.json"
    shell:
        """
        cat {input} > {output}
        """
