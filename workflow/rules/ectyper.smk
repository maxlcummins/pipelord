import re
import subprocess
import os

print(os.getcwd())

rule ectyper_run:
    input:
        assembly = config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta"
    output:
        directory(config['outdir']+"/{prefix}/ectyper/{sample}"),
        config['outdir']+"/{prefix}/ectyper/{sample}/output.tsv"
    conda:
        "../envs/ectyper.yaml"
    threads:
        1
    log:
        "logs/{prefix}/ectyper/{sample}_ectyper.log"
    shell:
        "ectyper -i {input} -o {output[0]} &> {log}"

rule run_ectyper_summarise:
    input:
        ectyper_summary=expand(config['outdir']+"/{prefix}/ectyper/{sample}/output.tsv", prefix=prefix, sample=sample_ids)
    output:
        combine_ectyper=config['outdir']+"/{prefix}/summaries/ectyper.txt"
    params:
        extra="",
    log:
        "logs/{prefix}/summaries/combine_ectyper.log",
    threads: 1
    script:
        "../../scripts/combine_ectyper.py"
