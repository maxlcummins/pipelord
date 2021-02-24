import re
import subprocess
import os

prefix = config['prefix']
maxthreads = snakemake.utils.available_cpu_count()

db_location = config['gene_db_location']
gene_dbs = expand(config['gene_dbs'])

logs = config['base_log_outdir']

rule mlst_run:
    input:
        assembly = config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta"
    output:
        config['outdir']+"/{prefix}/mlst/{sample}_mlst.txt"
    conda:
        "../envs/mlst.yaml"
    threads:
        1
    shell:
        "mlst {input} > {output}"
