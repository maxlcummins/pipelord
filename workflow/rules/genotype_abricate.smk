import re
import subprocess
import os

prefix = config['prefix']
maxthreads = snakemake.utils.available_cpu_count()

db_location = config['gene_db_location']
gene_dbs = expand(config['gene_dbs'])

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
