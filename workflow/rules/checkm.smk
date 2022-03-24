import re
import subprocess
import os

prefix = config["prefix"]
outdir = config["outdir"]
maxthreads = snakemake.utils.available_cpu_count()

logs = config["base_log_outdir"]

if platform.system() == "Darwin":
    checkm_env = "../envs/checkm_macOSX.yaml"
else:
    checkm_env = "../envs/checkm.yaml"

if platform.system() == "Darwin":
    rule install_pplacer:
        output:
            dummy_out = temp("downloaded_pplacer")
        conda:
            checkm_env
        threads:
            maxthreads
        shell:
            """
            FILE=$CONDA_PREFIX/bin/pplacer
            if test -f "$FILE"; then
                echo "$FILE exists. pplacer will not be downloaded"
            else
                echo "Downloading pplacer..."
                wget -q https://github.com/matsen/pplacer/releases/download/v1.1.alpha17/pplacer-Darwin-v1.1.alpha17.zip
                echo "Unzipping pplacer..."
                unzip pplacer-Darwin-v1.1.alpha17.zip && rm pplacer-Darwin-v1.1.alpha17.zip
                echo "Moving pplacer and associated scripts to $CONDA_PREFIX/bin"
                mv pplacer-Darwin-v1.1.alpha17*/* $CONDA_PREFIX/bin && rmdir pplacer-Darwin-v1.1.alpha17*/
                ls $CONDA_PREFIX/bin/pplacer
            fi
            touch {output}
            """

rule checkm_tree_and_tree_qa:
    input:
        dummy_out = "downloaded_pplacer",
        assemblies = config['outdir']+"/{prefix}/shovill/assemblies",
    output:
        directory(config['outdir']+"/{prefix}/QC_workflow/checkm")
    conda:
        checkm_env
    log:
        tree = config['base_log_outdir']+"/{prefix}/QC_workflow/checkm/checkm_tree_and_tree_qa.log"
    threads:
        maxthreads
    shell:
        """
        checkm tree {input.assemblies} -x fasta {output} -t {threads}
        checkm tree_qa {output}
        """

rule checkm_lineage_set:
    input:
        config['outdir']+"/{prefix}/QC_workflow/checkm"
    output:
        markers = config['outdir']+"/{prefix}/QC_workflow/checkm_markers/markers"
    conda:
        checkm_env
    log:
        config['base_log_outdir']+"/{prefix}/QC_workflow/checkm/checkm_lineage_set.log"
    threads:
        maxthreads
    shell:
        "checkm lineage_set {input} {output}"

rule checkm_analyze:
    input:
        assemblies = config['outdir']+"/{prefix}/shovill/assemblies",
        markers = config['outdir']+"/{prefix}/QC_workflow/checkm_markers/markers"
    output:
        temporary(config['outdir']+"/{prefix}/QC_workflow/checkm_dummy/alignment_info.tsv")
    conda:
        checkm_env
    log:
        config['base_log_outdir']+"/{prefix}/QC_workflow/checkm/checkm_analyze.log"
    threads:
        maxthreads
    shell:
        """
        checkm analyze {input.markers} {input.assemblies} {output} -t {threads} -x fasta
        touch {output}
        """

rule checkm_qa:
    input:
        markers = config['outdir']+"/{prefix}/QC_workflow/checkm_markers/markers",
        pseudoinput = config['outdir']+"/{prefix}/QC_workflow/checkm_dummy/alignment_info.tsv"
    output:
        config['outdir']+"/{prefix}/QC_workflow/checkm_qa/qa.tsv"
    conda:
        checkm_env
    log:
        config['base_log_outdir']+"/{prefix}/checkm/checkm_qa.log"
    shell:
        "checkm qa {input.checkmdir} {output} -f {output} -o 2 --tab_table"


ruleorder: install_pplacer > checkm_tree_and_tree_qa > checkm_lineage_set > checkm_analyze > checkm_qa