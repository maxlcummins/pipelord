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

#Download the appropriate blast tool depending on the OS
if path.exists("resources/tools/ncbi-blast-2.13.0+") == False:
    print('blast directory not located, downloading blast...')
    if platform.system() == "Linux":
        os.system("wget -O resources/tools/ncbi-blast-2.13.0+.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.12.0/ncbi-blast-2.12.0+-x64-linux.tar.gz")
    elif platform.system() == "Windows":
        os.system("wget -O resources/tools/ncbi-blast-2.13.0+.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.13.0+-x64-win64.tar.gz")
    elif platform.system() == "Darwin":
        os.system("wget -O resources/tools/ncbi-blast-2.13.0+.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.13.0+-x64-macosx.tar.gz")
    print('Unzipping...')
    os.system("tar -xvf resources/tools/ncbi-blast-2.13.0+.tar.gz -C resources/tools/")

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

if config["genotype_modules"]["run_genome_assembly"] == False:
    (sample_ids,) = glob_wildcards(config["genome_path"]+"/{sample}.fasta")
else:
    (sample_ids,) = glob_wildcards(config["genome_path"]+"/{sample}.R1.fastq.gz")

print(sample_ids,)

rule all:
    input:
        expand(config["outdir"]+"/{prefix}/fastp/{sample}.R1.fastq.gz",sample=sample_ids,prefix=prefix) if config["qc_modules"]["run_fastp"] else [],
        expand(config['outdir']+"/{prefix}/summaries/bracken_report.txt", prefix=prefix) if config["qc_modules"]["run_kraken2_and_bracken"] else [],
        expand(config['outdir']+"/{prefix}/QC_workflow/checkm_qa/qa.tsv", prefix=prefix) if config ["qc_modules"]["run_checkm"] and platform.system() == "Linux" else [],
        expand(config["outdir"]+"/{prefix}/QC_workflow/summaries/gunc_report.txt", prefix=prefix) if config ["qc_modules"]["run_gunc"] else[],
        expand(config["outdir"]+"/{prefix}/shovill/assembly_stats/{sample}_assembly_stats.txt", sample=sample_ids, prefix=prefix) if config["qc_modules"]["run_assembly_stats"] else [],
        expand( config["outdir"]+"/{prefix}/shovill/assemblies/{sample}.fasta", sample=sample_ids, prefix=prefix) if config["qc_modules"]["run_genome_assembly"] else [],


onsuccess:
#Announce job finish
    print("Job finished")
    #Email user
    if config['email'] != '':
        shell("echo Your Snakemake job with prefix \'{prefix}\' has finished. It has been written to \'{outdir}/{prefix}/summaries/\' | mail -s 'Snakemake job has finished' {email}")

if config["qc_modules"]["run_checkm"]:
    include: "rules/checkm.smk"
if config["qc_modules"]["run_kraken2_and_bracken"]:
    include: "rules/species_Id.smk"
if config["qc_modules"]["run_gunc"]:
    include: "rules/gunc.smk"
if config["qc_modules"]["run_fastp"]:
    include: "rules/read_cleaning.smk"
if config["qc_modules"]["run_assembly_stats"] or config["qc_modules"]["run_genome_assembly"]:
    include: "rules/genome_assembly.smk"