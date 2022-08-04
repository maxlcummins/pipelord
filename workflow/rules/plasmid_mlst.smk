import re
import subprocess
import os

if path.exists("resources/tools/pmlst") == False:
    print('pmlst directory not located, downloading pmlst...')
    os.system("git clone https://git@bitbucket.org/genomicepidemiology/pmlst.git resources/tools/pmlst")
    os.system("git clone https://git@bitbucket.org/genomicepidemiology/pmlst_db.git resources/tools/pmlst/pmlst_db")


prefix = config['prefix']
scheme = config['pmlst_scheme']


logs = config['base_log_outdir']

rule pMLST_run:
    input:
        assembly = config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta"
    output:
        results = config['outdir']+"/{prefix}/pMLST/{scheme}/{sample}.out/results.txt",
    conda:
        "../envs/pMLST.yaml"
    threads: 4
    params:
        output_dir = config['outdir']+"/{prefix}/pMLST/{scheme}/{sample}.out",
        tmp = config['outdir']+"/{prefix}/pMLST/{scheme}/temp_{sample}",
        pmlst_tool = config['pmlst_tool'],
        pmlst_script_path = config['pMLST_script_path'],
        pmlst_db_path = config['pMLST_db_path'],
        db = "{scheme}",
        scheme = scheme,
    shell:
        """
        python3 {params.pmlst_script_path} -i {input} -o {params.output_dir} -p {params.pmlst_db_path} {params.pmlst_tool} $CONDA_PREFIX/bin/blastn -s {params.db} -x -t {params.tmp}
        rm -rf {params.tmp}
        """

rule run_pmlst_summarise:
    input:
        pMLST_summary=expand(config['outdir']+"/{prefix}/pMLST/{scheme}/{sample}.out/results.txt", prefix=prefix, scheme=scheme, sample=sample_ids)
    output:
        combine_pMLST=config["outdir"]+"/{prefix}/summaries/pMLST.txt",
    params:
        extra="",
    log:
        "logs/{prefix}/summaries/combine_pMLST.log",
    conda:
        "../envs/abritamr.yaml"
    threads: 1
    script:
        "../../scripts/combine_pMLST.py"

if 'incf' in config['pmlst_scheme']:
    rule incf_touch:
        input:
            incf_results = config['outdir']+"/{prefix}/pMLST/incf/{sample}.out/results.txt",
        output:
            incf_tab = config['outdir']+"/{prefix}/pMLST/incf/{sample}.out/results_tab.tsv"
        threads: 1
        shell:
            """
            touch {input}
            """

if 'incf' in config['pmlst_scheme']:
    rule run_incf_summarise:
        input:
            pMLST_summary=expand(config['outdir']+"/{prefix}/pMLST/{scheme}/{sample}.out/results.txt", prefix=prefix, scheme=scheme, sample=sample_ids)
        output:
            IncF_RST_simple=config["outdir"]+"/{prefix}/summaries/IncF_RST.txt",
            IncF_RST_full=config["outdir"]+"/{prefix}/summaries/IncF_RST_full.txt"
        params:
            extra="",
        log:
            "logs/{prefix}/summaries/combine_incF.log",
        conda:
            "../envs/abritamr.yaml"
        threads: 1
        script:
            "../../scripts/combine_IncF_RST.py"