import re
import subprocess
import os

prefix = config["prefix"]
outdir = config["outdir"]
maxthreads = snakemake.utils.available_cpu_count()

rule QC_summarise:
    input:
        bracken_report = config['outdir']+"/{prefix}/QC_workflow/summaries/bracken_report.txt" if config["qc_modules"]["run_kraken2_and_bracken"] else [],
        checkm_report = config['outdir']+"/{prefix}/QC_workflow/summaries/checkm_qa.txt" if config ["qc_modules"]["run_checkm"] and platform.system() == "Linux" else [],
        gunc_report = config["outdir"]+"/{prefix}/QC_workflow/summaries/gunc_report.txt" if config ["qc_modules"]["run_gunc"] else[],
        assembly_stats_report = config["outdir"]+"/{prefix}/QC_workflow/summaries/assembly_stats.txt"
    output:
        config['outdir']+"/{prefix}/QC_workflow/summaries/QC_report.txt"
    log:
        "logs/{prefix}/summaries/summarise_qc.log",
    params:
        names = sample_ids,
    script:
        "scripts/summarise_qc.py"