import re
import subprocess
import os
from os import path
import git
import re

configfile: "config/config_template.yaml"

outdir = config["outdir"]
prefix = config["prefix"]
gene_dbs = expand(config["gene_dbs"])
plasmid_screen_db = expand(config["plasmid_screen_db"])

scheme = config["pmlst_scheme"]
email = config["email"]

if path.exists("resources/tools") == False:
    print("tools directory not located, creating tools directory...")
    os.system("mkdir -p resources/tools")

if config["input_type"] == "assemblies":
    (sample_ids,) = glob_wildcards(config["genome_path"]+"/{sample}.fasta")
elif config["input_type"] == "reads":
    (sample_ids,) = glob_wildcards(config["genome_path"]+"/{sample}.R1.fastq.gz")
else: "Config variable 'input_type' must be either 'assemblies' or 'reads'. Please check the config file"

if config["input_type"] == "assemblies":
    (sample_ids,) = glob_wildcards(config["genome_path"]+"/{sample}.fasta")
    print("Detected "+str(len(sample_ids))+" genomes")
    for i in sample_ids:
        if os.path.getsize(config["genome_path"]+"/"+i+".fasta") < 1000000:
            sample_ids.remove(i)
            print(i+" fasta too small, excluded from analysis")
    print(str(len(sample_ids))+" genomes of an appropriate size for inclusion... (assembly over 1MB)")
elif config["input_type"] == "reads":
    (sample_ids,) = glob_wildcards(config["genome_path"]+"/{sample}.R1.fastq.gz")
    print("Detected "+str(len(sample_ids))+" genomes")
    for i in sample_ids:
        if os.path.getsize(config["genome_path"]+"/"+i+".R1.fastq.gz") < 20000000:
            print(i+" R1 too small, excluded from analysis")
            sample_ids.remove(i)
        elif os.path.getsize(config["genome_path"]+"/"+i+".R2.fastq.gz") < 20000000:
            print(i+" R2 too small, excluded from analysis")
            sample_ids.remove(i)
    print(str(len(sample_ids))+" genomes of an appropriate size for inclusion... (R1 and R2 over 20MB each)")
else: print("Config variable 'input_type' must be either 'assemblies' or 'reads'. Please check the config file")

print(sample_ids,)

#Access our files of file names for our subsets
subset_fofn = config["subset_fofn"]

#Access our names for our subset groups
subset_prefixes = config["subset_name"]

rule all:
    input: expand(config['outdir']+"/{prefix}/summaries/{prefix}.distances.tsv", prefix=prefix),
    #input: expand(config['outdir']+"/{prefix}/summaries/{prefix}_100_SNPs.distances.tsv", prefix=prefix),

if config['input_type'] == "reads":
    rule ska_fastq:
        input:
            r1_filt = config['outdir']+"/{prefix}/fastp/{sample}.R1.fastq.gz",
            r2_filt = config['outdir']+"/{prefix}/fastp/{sample}.R2.fastq.gz"
        output:
            config['outdir']+"/{prefix}/ska/fastq/{sample}.skf"
        conda:
            "../../envs/ska.yaml"
        log:
            config['base_log_outdir']+"/{prefix}/ska/fastq/{sample}.log"
        params:
            out = config['outdir']+"/{prefix}/ska/fastq/{sample}",
            log_dir = config['base_log_outdir']+"/{prefix}/ska/fastq"
        threads: 1
        shell:
            """
            mkdir -p {params.log_dir}
            ska fastq \
                -o {params.out} \
                {input.r1_filt} \
                {input.r2_filt} 2> {log}
            """

if config['input_type'] == "assemblies":
    rule ska_fasta:
        input:
            assembly = config['genome_path']+"/{sample}.fasta"
        output:
            config['outdir']+"/{prefix}/ska/fasta/{sample}.skf"
        conda:
            "../../envs/ska.yaml"
        log:
            config['base_log_outdir']+"/{prefix}/ska/fasta/{sample}.log"
        params:
            out = config['outdir']+"/{prefix}/ska/fasta/{sample}",
            log_dir = config['base_log_outdir']+"/{prefix}/ska/fasta"
        threads: 1
        shell:
            """
            mkdir -p {params.log_dir}
            ska fasta \
                -o {params.out} \
                {input.assembly} 2> {log}
            """

rule ska_distance:
    input:
        expand(config['outdir']+"/{prefix}/ska/fasta/{sample}.skf", prefix=prefix, sample=sample_ids) if config['input_type'] == "assemblies" else [],
        expand(config['outdir']+"/{prefix}/ska/fastq/{sample}.skf", prefix=prefix, sample=sample_ids) if config['input_type'] == "reads" else []
    output:
        config['outdir']+"/{prefix}/summaries/{prefix}.clusters.tsv",
        config['outdir']+"/{prefix}/summaries/{prefix}.distances.tsv"
    conda:
        "../../envs/ska.yaml"
    log:
        config['base_log_outdir']+"/{prefix}/ska/distance.log"
    params:
        out = expand(config['outdir']+"/{prefix}/summaries/{prefix}", prefix=prefix, sample=sample_ids),
    threads: 30
    shell:
        """
        ska distance \
            -s 25 \
            -i 0.95 \
            -o {params.out} \
            {input}  2> {log}
        """

rule ska_distance_far:
    input:
        expand(config['outdir']+"/{prefix}/ska/fasta/{sample}.skf", prefix=prefix, sample=sample_ids) if config['input_type'] == "assemblies" else [],
        expand(config['outdir']+"/{prefix}/ska/fastq/{sample}.skf", prefix=prefix, sample=sample_ids) if config['input_type'] == "reads" else []
    output:
        config['outdir']+"/{prefix}/summaries/{prefix}_100_SNPs.clusters.tsv",
        config['outdir']+"/{prefix}/summaries/{prefix}_100_SNPs.distances.tsv"
    conda:
        "../../envs/ska.yaml"
    log:
        config['base_log_outdir']+"/{prefix}/ska/distance_100_SNPs.log"
    params:
        out = expand(config['outdir']+"/{prefix}/summaries/{prefix}_100_SNPs", prefix=prefix, sample=sample_ids),
    threads: 30
    shell:
        """
        ska distance \
            -s 100 \
            -i 0.95 \
            -o {params.out} \
            {input}  2> {log}
        """

onerror:
    print("An error occurred")
    # Email user
    if config["email"] != "":
        shell(
            "echo Your Snakemake SKA job with prefix '{prefix}' has failed. | mail -s 'Snakemake SKA job has failed' {email}"
        )

