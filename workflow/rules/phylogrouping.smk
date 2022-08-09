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

if path.exists("resources/tools/clermontyping") == False:
    print('clermontyping tool not located, downloading clermontyping...')
    rule clermontyping_download:
        output:
            directory("resources/tools/clermontyping")
        conda:
            "../envs/clermontyping.yaml"
        threads:
            1
        shell:
            """
            git clone https://github.com/A-BN/ClermonTyping.git resources/tools/clermontyping
            """

rule clermontyping_run:
    input:
        github_repo = "resources/tools/clermontyping",
        assembly = config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta"
    output:
        config['outdir']+"/{prefix}/clermontyping/{sample}/{sample}_phylogroups.txt",
    conda:
        "../envs/clermontyping.yaml"
    threads:
        1
    log:
        config['base_log_outdir']+"/{prefix}/clermontyping/{sample}.log"
    params:
        config['outdir']+"/{prefix}/clermontyping/{sample}"
    shell:
        """
        sh {input.github_repo}/clermonTyping.sh --fasta {input.assembly} --name {wildcards.sample} --minimal &> {log}
        ls {params}
        mv {wildcards.sample}/* {params}
        rmdir {wildcards.sample}
        """

rule run_clermontyping_summarise:
    input:
        clermontyping = expand(config['outdir']+"/{prefix}/clermontyping/{sample}/{sample}_phylogroups.txt", prefix=prefix, sample=sample_ids)
    output:
        config['outdir']+"/{prefix}/summaries/clermontyping.txt"
    log:
        "logs/{prefix}/summaries/combine_clermontyping.log",
    conda:
        "../envs/abritamr.yaml" #Need to replace with one with pandas
    threads: 1
    script:
        "../../scripts/combine_clermontyping.py"