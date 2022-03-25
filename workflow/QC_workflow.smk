import re
import subprocess
import os
from os import path
import git
import platform


configfile: "config/config_template.yaml"

if platform.system() == "Darwin" and config ["qc_modules"]["run_checkm"]:
    print("CheckM enabled but operating system detected as Darwin (MacOSX). CheckM is not supported on Mac, therefore it will not be run")
    print("If you want to run CheckM you must run this pipeline on a Linux system. A future release will support a containerised, Linux-based CheckM enironment.")

if path.exists("resources/tools/kma") == False:
    print('kma directory not located, downloading kma...')
    os.system("git clone https://bitbucket.org/genomicepidemiology/kma.git resources/tools/kma")
    os.system("cd resources/tools/kma && make")
    os.system("cd ../../..")

outdir = config["outdir"]
prefix = config["prefix"]

logs = config["base_log_outdir"]

email = config["email"]

if path.exists("resources/tools") == False:
    print("tools directory not located, creating tools directory...")
    os.system("mkdir -p resources/tools")

if config["input_type"] == "assemblies":
    (sample_ids,) = glob_wildcards(config["genome_path"]+"/{sample}.fasta")
elif config["input_type"] == "reads":
    (sample_ids,) = glob_wildcards(config["genome_path"]+"/{sample}.R1.fastq.gz")
else: "Config variable 'input_type' must be either 'assemblies' or 'reads'. Please check the config file"

#print("\nFrom config file, input type selected as '"+config["input_type"]+"'\n")
#print("Genomes detected:")
#print(sample_ids,)
#print("\n")

rule all:
    input:
        expand(config['outdir']+"/{prefix}/QC_workflow/summaries/multiqc_fastp.txt", prefix=prefix) if config["qc_modules"]["run_multiqc_fastp"] else [],
        expand(config['outdir']+"/{prefix}/QC_workflow/summaries/QC_report.txt",prefix=prefix) if config["qc_modules"]["run_qc_summary"] else [],
        expand(config["outdir"]+"/{prefix}/fastp/{sample}.R1.fastq.gz",sample=sample_ids,prefix=prefix) if config["input_type"] == "reads" else [],
        expand(config['outdir']+"/{prefix}/QC_workflow/summaries/bracken_report.txt", prefix=prefix) if config["qc_modules"]["run_kraken2_and_bracken"] else [],
        expand(config['outdir']+"/{prefix}/QC_workflow/bracken/{sample}.bracken.txt", sample=sample_ids, prefix=prefix) if config["qc_modules"]["run_kraken2_and_bracken"] else [],
        expand(config['outdir']+"/{prefix}/QC_workflow/summaries/checkm_qa.txt", prefix=prefix) if config ["qc_modules"]["run_checkm"] and platform.system() == "Linux" else [],
        expand(config["outdir"]+"/{prefix}/QC_workflow/summaries/gunc_report.txt", prefix=prefix) if config ["qc_modules"]["run_gunc"] else[],
        expand(config["outdir"]+"/{prefix}/QC_workflow/summaries/assembly_stats.txt",prefix=prefix) if config["qc_modules"]["run_assembly_stats"] else [],
        expand( config["outdir"]+"/{prefix}/shovill/assemblies/{sample}.fasta", sample=sample_ids, prefix=prefix) if config["qc_modules"]["run_genome_assembly"] else [],


onsuccess:
#Announce job finish
    print("Job finished")
    #Email user
    if config['email'] != '':
        shell("echo Your Snakemake job with prefix \'{prefix}\' has finished. It has been written to \'{outdir}/{prefix}/summaries/\' | mail -s 'Snakemake job has finished' {email}")

include: "rules/read_cleaning.smk"

if config["qc_modules"]["run_checkm"]:
    include: "rules/checkm.smk"
if config["qc_modules"]["run_kraken2_and_bracken"]:
    include: "rules/species_id.smk"
if config["qc_modules"]["run_gunc"]:
    include: "rules/gunc.smk"
if config["qc_modules"]["run_assembly_stats"] or config["qc_modules"]["run_genome_assembly"]:
    include: "rules/genome_assembly.smk"
if config["qc_modules"]["run_qc_summary"]:
    include: "rules/QC_summary.smk"
if config["qc_modules"]["run_multiqc_fastp"]:
    include: "rules/multiqc_fastp.smk"