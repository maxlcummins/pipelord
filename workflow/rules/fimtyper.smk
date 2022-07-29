import re
import subprocess
import os
from os import path
import git

prefix = config['prefix']

if path.exists("resources/tools/fimtyper") == False:
    print('fimtyper directory not located, downloading fimtyper...')
    os.system("git clone https://git@bitbucket.org/genomicepidemiology/fimtyper.git resources/tools/fimtyper")
    os.system("git clone https://git@bitbucket.org/genomicepidemiology/fimtyper_db.git resources/tools/fimtyper/fimtyper_db")
    #os.system("perl -p -i -e 's@^sample_name = filename.split.*@sample_name = filename.split(\".\")[0]@g' resources/tools/fimtyper/fimtyper.py")

rule validate_db:
    input:
        "resources/tools/fimtyper/fimtyper_db"
    output:
        config['outdir']+"/{prefix}/fimtyper/dummy_out.txt"
    conda:
        "../envs/fimtyper.yaml"
    log:
        config['base_log_outdir']+"/{prefix}/fimtyper/fimtyper_db_check.log"
    threads: 4
    params:

    shell:
        """
        python3 resources/tools/fimtyper/VALIDATE_DB {input} 2> {log}
        touch {output}
        """

rule fimtyper_run:
    input:
        dummy = config['outdir']+"/{prefix}/fimtyper/dummy_out.txt",
        assembly = config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta"
    output:
        raw = config['outdir']+"/{prefix}/fimtyper/{sample}/results_tab.txt",
        renamed = config['outdir']+"/{prefix}/fimtyper/{sample}/{sample}.txt"
    conda:
        "../envs/fimtyper.yaml"
    log:
        config['base_log_outdir']+"/{prefix}/fimtyper/fimtyper_run/{sample}.log"
    threads: 4
    params:
        output_dir = config['outdir']+"/{prefix}/fimtyper/{sample}",
        blastn_path = config['blast_bin']
    shell:
        """
        perl resources/tools/fimtyper/fimtyper.pl -d resources/tools/fimtyper/fimtyper_db -b $CONDA_PREFIX -i {input.assembly} -o {params.output_dir} -k 95.0 -l 0.6 &> {log}
        cp {output.raw} {output.renamed}
        """

rule name_append_fimtyper:
    input:
        config['outdir']+"/{prefix}/fimtyper/{sample}/{sample}.txt"
    output:
        trimmed = config['outdir']+"/{prefix}/fimtyper/{sample}/{sample}_trimmed.txt",
        named = config['outdir']+"/{prefix}/fimtyper/{sample}/{sample}_named.txt"
    threads: 1
    shell:
        """
        sed '0,/^FimH type -/d' {input} > {output.trimmed}
        awk 'NR == 1 {{print "name\t" $0; next;}}{{print FILENAME "\t" $0;}}' {output.trimmed} > {output.named}
        """

rule run_fimtyper_summarise:
    input:
        fimtyper_summary=expand(config['outdir']+"/{prefix}/fimtyper/{sample}/{sample}_named.txt", prefix=prefix, sample=sample_ids)
    output:
        combine_fimtyper=config["outdir"]+"/{prefix}/summaries/fimtyper.txt"
    params:
        extra="",
    log:
        "logs/{prefix}/summaries/combine_fimtyper.log",
    threads: 1
    script:
        "../../scripts/combine_fimtyper.py"