import re
import subprocess
import os

prefix = config['prefix']

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

rule run_mlst_summarise:
    input:
        mlst_summary=expand(config['outdir']+"/{prefix}/mlst/{sample}_mlst.txt", prefix=prefix, sample=sample_ids)
    output:
        combine_mlst=config['outdir']+"/{prefix}/summaries/mlst.txt"
    params:
        extra="",
    log:
        "logs/{prefix}/summaries/combine_mlst.log",
    threads: 1
    script:
        "../../scripts/combine_mlst.py"
