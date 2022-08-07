import re
import subprocess
import os

rule ezclermont_run:
    input:
        assembly = config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta"
    output:
        simple = config['outdir']+"/{prefix}/ezclermont/{sample}/ezclermont.txt",
        full = config['outdir']+"/{prefix}/ezclermont/{sample}/ezclermont.log"
    conda:
        "../envs/ezclermont.yaml"
    threads:
        1
    params:
        ezclermont_outdir = config["outdir"]+"/{prefix}/ezclermont/{sample}"
    shell:
        """
        mkdir -p {params.ezclermont_outdir}
        ! ezclermont {input} --logfile {output.full} > {output.simple}
        """

rule run_ezclermont_summarise:
    input:
        ezclermont_summary=expand(config['outdir']+"/{prefix}/ezclermont/{sample}/ezclermont.log", prefix=prefix, sample=sample_ids)
    output:
        combine_ezclermont=config['outdir']+"/{prefix}/summaries/ezclermont.txt"
    log:
        "logs/{prefix}/summaries/combine_ezclermont.log",
    conda:
        "../envs/abritamr.yaml" #Need to replace with one with pandas
    threads: 1
    script:
        "../../scripts/combine_ezclermont.py"
