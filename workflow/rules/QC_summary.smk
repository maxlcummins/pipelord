import re
import subprocess
import os

prefix = config["prefix"]
outdir = config["outdir"]
maxthreads = snakemake.utils.available_cpu_count()

rule QC_summarise:
    input:
        bracken_report = 
        assembly_stats_repor = 
        gunc_report = 
        checkm_report = 
    output:
        config['outdir']+"/{prefix}/QC_workflow/summaries/QC_report.txt"
    conda:
        "../envs/gunc.yaml"
    script:
        "scripts/summarise_qc.py""