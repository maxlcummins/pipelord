import re
import subprocess
import os

prefix = config["prefix"]
outdir = config["outdir"]
maxthreads = snakemake.utils.available_cpu_count()

logs = config["base_log_outdir"]

rule multiqc:
    input:
        json = expand(config['outdir']+"/{prefix}/fastp/{sample}.fastp.json", prefix=prefix, sample=sample_ids)
    output:
        multiqc_out = directory(config['outdir']+"/{prefix}/QC_workflow/multiqc"),
        summary = config['outdir']+"/{prefix}/QC_workflow/summaries/multiqc_fastp.txt"
    conda:
        "../envs/multiqc.yaml"
    log:
        config['base_log_outdir']+"/{prefix}/QC_workflow/multiqc/multiqc.log"
    threads:
        1
    shell:
        """
        mutliqc {input} -o {output}
        cp {output}/multiqc_general_stats.txt {output.summary}
        """
