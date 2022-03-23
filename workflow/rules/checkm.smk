import re
import subprocess
import os

prefix = config["prefix"]
outdir = config["outdir"]
maxthreads = snakemake.utils.available_cpu_count()

logs = config["base_log_outdir"]

rule checkm_tree_and_tree_qa:
    input:
        config['outdir']+"/{prefix}/shovill/assemblies"
    output:
        directory(config['outdir']+"/{prefix}/checkm")
    conda:
        "../envs/checkm.yaml"
    log:
        tree = config['base_log_outdir']+"/{prefix}/checkm/checkm_tree_and_tree_qa.log"
    threads:
        maxthreads
    shell:
        """
        checkm tree {input} -x fasta {output} -t {threads}
        checkm tree_qa {output}
        """

rule checkm_lineage_set:
    input:
        config['outdir']+"/{prefix}/checkm"
    conda:
        "../envs/checkm.yaml"
    log:
        config['base_log_outdir']+"/{prefix}/checkm/checkm_lineage_set.log"
    threads:
        maxthreads
    shell:
        "checkm lineage_set {input} {input}/markers"

rule checkm_analyze:
    input:
        assemblies = config['outdir']+"/{prefix}/shovill/assemblies",
        checkmdir = config['outdir']+"/{prefix}/checkm"
    output:
        directory(config['outdir']+"/{prefix}/checkm")
    conda:
        "../envs/checkm.yaml"
    log:
        config['base_log_outdir']+"/{prefix}/checkm/checkm_analyze.log"
    threads:
        maxthreads
    shell:
        "checkm analyze {input.checkmdir}/markers {input.assemblies} {output} -t {threads} -x fasta"

rule checkm_qa:
    input:
        checkmdir = config['outdir']+"/{prefix}/checkm"
    output:
        config['outdir']+"/{prefix}/checkm_qa/qa.tsv"
    conda:
        "../envs/checkm.yaml"
    log:
        config['base_log_outdir']+"/{prefix}/checkm/checkm_qa.log"
    shell:
        "checkm qa {input.checkmdir}/markers {output} -f {output} -o 2 --tab_table"


ruleorder: checkm_tree_and_tree_qa > checkm_lineage_set > checkm_analyze > checkm_qa