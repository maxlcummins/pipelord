import re
import subprocess
import os

prefix = config['prefix']

abritamr_outdir = config["outdir"]+"/"+config['prefix']+"/abritamr",

logs = config['base_log_outdir']

rule update_database:
    output:
        temp(config['outdir']+"/{prefix}/abritamr/dummy_out.txt")
    conda:
        "../envs/abritamr.yaml"
    threads:
        1
    log:
        config['base_log_outdir']+"/{prefix}/abritamr/update_database.log"
    shell:
        """
        amrfinder -U &> {log}
        touch {output}
        """

rule run_abritamr:
    input:
        contigs = config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta",
        dummy_db_file = config['outdir']+"/{prefix}/abritamr/dummy_out.txt"
    output:
        resistance = config["outdir"]+"/{prefix}/abritamr/{sample}/summary_matches.txt",
        partials = config["outdir"]+"/{prefix}/abritamr/{sample}/summary_partials.txt",
        virulence = config["outdir"]+"/{prefix}/abritamr/{sample}/summary_virulence.txt",
        amrfinder = config["outdir"]+"/{prefix}/abritamr/{sample}/amrfinder.out"
    params:
        abritamr_out = abritamr_outdir,
        species = config['abritamr_species']
    log:
        config['base_log_outdir']+"/{prefix}/abritamr/run/{sample}.log"
    threads: 4
    conda:
        "../envs/abritamr.yaml"
    shell:
        """
        abritamr run --prefix {params.abritamr_out}/{wildcards.sample} --contigs {input.contigs}  --jobs {threads} --species {params.species} &> {log}       
        """

rule run_abritamr_summarise:
    input:
        abritamr_resistance=expand(config["outdir"]+"/{prefix}/abritamr/{sample}/summary_matches.txt", prefix=prefix, sample=sample_ids),
        abritamr_virulence=expand(config["outdir"]+"/{prefix}/abritamr/{sample}/summary_virulence.txt", prefix=prefix, sample=sample_ids),
        abritamr_partials=expand(config["outdir"]+"/{prefix}/abritamr/{sample}/summary_partials.txt", prefix=prefix, sample=sample_ids)
    output:
        abritamr_resistance=config["outdir"]+"/{prefix}/summaries/abritamr_resistance.txt",
        abritamr_virulence=config["outdir"]+"/{prefix}/summaries/abritamr_virulence.txt",
        abritamr_partials=config["outdir"]+"/{prefix}/summaries/abritamr_partials.txt"
    log:
        config['base_log_outdir']+"/{prefix}/abritamr/summarise.log"
    threads: 1
    script:
        "../../scripts/combine_abritamr.py"

rule run_abritamr_raw_amrfinder_summarise:
    input:
        amrfinder=expand(config["outdir"]+"/{prefix}/abritamr/{sample}/amrfinder.out", prefix=prefix, sample=sample_ids)
    output:
        abritamr_resistance=config["outdir"]+"/{prefix}/summaries/amrfinder_raw.txt"
    log:
        "logs/{prefix}/summaries/combine_abritamr.log",
    threads: 1
    script:
        "../../scripts/combine_abritamr_amrfinder_raw.py"

ruleorder: update_database > run_abritamr