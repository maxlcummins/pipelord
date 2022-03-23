import re
import subprocess
import os

prefix = config['prefix']
maxthreads = snakemake.utils.available_cpu_count()

db_location = config['gene_db_location']
gene_dbs = expand(config['gene_dbs'])
plasmid_screen_db = expand(config['plasmid_screen_db'])

logs = config['base_log_outdir']


rule abricate_run:
    input:
        assembly = config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta"
    output:
        tab = config['outdir']+"/{prefix}/abricate/{gene_db}/{sample}.tab",
    log:
        config['base_log_outdir']+"/{prefix}/abricate/run/{gene_db}/{sample}.log"
    conda:
        "../envs/abricate.yaml"
    threads:
        6
    params:
        db = "{gene_db}",
        gene_db = gene_dbs,
        datadir = config['gene_db_location']
    shell:
        "abricate --nopath --datadir {params.datadir} --db {params.db} {input.assembly} > {output} 2> {log}"

rule abricate_plasmid_run:
    input:
        assembly = config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta"
    output:
        config['outdir']+"/{prefix}/abricate_plasmids/{plasmid_screen_db}/{sample}.tab",
    log:
        config['base_log_outdir']+"/{prefix}/abricate/run/{plasmid_screen_db}/{sample}.log"
    conda:
        "../envs/abricate_plasmid.yaml"
    threads:
        6
    params:
        db = "{plasmid_screen_db}",
        plasmid_screen_db = plasmid_screen_db,
        datadir = config['gene_db_location']
    shell:
        "abricate --nopath --datadir {params.datadir} --db {params.db} {input.assembly} > {output} 2> {log}"

rule abricate_plasmid_combine:
    input:
        expand(config['outdir']+"/{prefix}/abricate/{plasmid_screen_db}/{sample}.tab", plasmid_screen_db=plasmid_screen_db, sample=sample_ids, prefix=prefix)
    output:
        out1 = temp(config['outdir']+"/{prefix}/abricate/plasmid_screen_tmp.txt"),
        out2 = temp(config['outdir']+"/{prefix}/abricate/plasmid_screen_tmp2.txt"),
        out = config['outdir']+"/{prefix}/abricate/plasmid_screen.txt",
    log:
        config['base_log_outdir']+"/{prefix}/abricate/run/plasmid_screen.log"
    conda:
        "abricate_plasmid.yaml"
    threads:
        1
    shell:
        """"
        cat {input} > {output.out1}
        perl -p -i -e 's@\.fasta@@g' {output.out1}
        #Rename first column name
        sed -E '1s/#FILE/name/g' {output.out1} > {output.out2}
        #Remove duplicate headers by negative grep
        grep -v '#FILE' {output.out2} > {output.out}
        """

rule run_abricate_summarise:
    input:
        abricate_summary=expand(config['outdir']+"/{prefix}/abricate/{gene_db}/{sample}.tab", prefix=prefix, gene_db=gene_dbs, sample=sample_ids)
    output:
        combine_abricate=config["outdir"]+"/{prefix}/summaries/genotype.txt"
    params:
        extra="",
    log:
        "logs/{prefix}/summaries/combine_abricate.log",
    threads: 1
    script:
        "../../scripts/combine_abricate.py"