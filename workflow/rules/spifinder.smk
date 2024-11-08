import re
import subprocess
import os
from os import path
import git

prefix = config['prefix']

if path.exists("resources/tools/spifinder") == False:
    print('spifinder directory not located, downloading spifinder...')
    os.system("git clone https://bitbucket.org/genomicepidemiology/spifinder.git resources/tools/spifinder")

if path.exists("/projects/AusGEM/databases/spifinder_db") == False:
    if path.exists("resources/tools/spifinder_db") == False:
        print('spifinder_db directory not located, downloading database...')
        os.system("git clone https://git@bitbucket.org/genomicepidemiology/spifinder_db.git resources/tools/spifinder_db")
        os.system("cd resources/tools/spifinder_db")
        os.system("python3 INSTALL.py kma_index")
        os.system("cd ../../..")


prefix = config['prefix']

logs = config['base_log_outdir']

if config['input_type'] == 'raw_reads':
    rule run_spifinder:
        input:
            r1_filt = config['outdir']+"/{prefix}/fastp/{sample}.R1.fastq.gz",
            r2_filt = config['outdir']+"/{prefix}/fastp/{sample}.R2.fastq.gz"
        output:
            spifinder_out = directory(config['outdir']+"/{prefix}/spifinder/spifinder_out/{sample}.out"),
            spifinder_summary = config['outdir']+"/{prefix}/spifinder/spi_summaries/{sample}.txt"
        log:
            out = config['base_log_outdir']+"/{prefix}/spifinder/run/{sample}_out.log",
            err = config['base_log_outdir']+"/{prefix}/spifinder/run/{sample}_err.log"
        threads:
            6
        conda:
            "../envs/spifinder.yaml"
        shell:
            """
            mkdir -p {output.spifinder_out}
            resources/tools/spifinder/spifinder.py -p resources/tools/spifinder_db -x -o {output.spifinder_out} -l 0.60 -t 0.95 -mp $CONDA_PREFIX/bin/kma -i {input.r1_filt} {input.r2_filt}
            cp {output.spifinder_out}/results_tab.tsv {output.spifinder_summary}
            """
    rule run_spifinder_summarise:
        input:
            #assembly-stats text files
            spifinder_summary=expand(config['outdir']+"/{prefix}/spifinder/spi_summaries/{sample}.txt", prefix=prefix, sample=sample_ids)
        output:
            #combined, tab-delimited assembly statistics
            combined_spifinder=config['outdir']+"/{prefix}/summaries/spi_report.txt"
        params:
            extra="",
        log:
            "logs/{prefix}/summaries/combine_spifinder.log",
        threads: 1
        script:
            "../../scripts/combine_spifinder.py"


if config['input_type'] == 'assemblies':
    rule run_spifinder:
        input:
            config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta"
        output:
            spifinder_out = directory(config['outdir']+"/{prefix}/spifinder/spifinder_out/{sample}.out"),
            spifinder_summary = config['outdir']+"/{prefix}/spifinder/spi_summaries/{sample}.txt"
        log:
            out = config['base_log_outdir']+"/{prefix}/spifinder/run/{sample}_out.log",
            err = config['base_log_outdir']+"/{prefix}/spifinder/run/{sample}_err.log"
        threads:
            5
        conda:
            "../envs/spifinder.yaml"
        shell:
            """
            mkdir -p {output.spifinder_out}
            resources/tools/spifinder/spifinder.py -p resources/tools/spifinder_db -x -o {output.spifinder_out} -l 0.60 -t 0.95 -mp $CONDA_PREFIX/bin/blastn -i {input}
            cp {output.spifinder_out}/results_tab.tsv {output.spifinder_summary}
            """
    rule run_spifinder_summarise:
        input:
            #assembly-stats text files
            spifinder_summary=expand(config['outdir']+"/{prefix}/spifinder/spi_summaries/{sample}.txt", prefix=prefix, sample=sample_ids)
        output:
            #combined, tab-delimited assembly statistics
            combine_spifinder=config['outdir']+"/{prefix}/summaries/spi_report.txt"
        params:
            extra="",
        log:
            "logs/{prefix}/summaries/combine_spifinder.log",
        threads: 1
        script:
            "../../scripts/combine_spifinder.py"

