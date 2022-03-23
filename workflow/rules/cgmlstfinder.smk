import re
import subprocess
import os
from os import path
import git

prefix = config['prefix']
maxthreads = snakemake.utils.available_cpu_count()

if config["cgmlst_schemes"] != "":
    if path.exists("resources/tools/cgmlstfinder") == False:
        print('cgmlstfinder directory not located, downloading cgmlstfinder...')
        os.system("git clone https://git@bitbucket.org/genomicepidemiology/cgmlstfinder.git resources/tools/cgmlstfinder")

if config["cgmlst_schemes"] != "":
    if path.exists("/projects/AusGEM/databases/cgmlstfinder_db/ecoli") == False:
        if path.exists("resources/dbs/cgmlstfinder_db/ecoli/ecoli.length.b") == False:
            print('cgmlstfinder_db directory not located, downloading database... This will take a long time but will only need to be done once')
            os.system("git clone https://git@bitbucket.org/genomicepidemiology/cgmlstfinder_db.git resources/dbs/cgmlstfinder_db")
            os.system("python resources/dbs/cgmlstfinder_db/INSTALL.py -s ecoli")

prefix = config['prefix']

logs = config['base_log_outdir']

if config['input_type'] == 'raw_reads':
    rule run_cgmlstfinder_assemblies:
        input:
            r1_filt = config['outdir']+"/{prefix}/fastp/{sample}.R1.fastq.gz",
            r2_filt = config['outdir']+"/{prefix}/fastp/{sample}.R2.fastq.gz"
        output:
            cgmlstfinder_out = directory(config['outdir']+"/{prefix}/cgmlstfinder/cgmlstfinder_out/{sample}.out"),
            cgmlst_summary = config['outdir']+"/{prefix}/cgmlstfinder/cgmlst_summaries/{sample}.txt"
        log:
            out = config['base_log_outdir']+"/{prefix}/cgmlstfinder/run/{sample}_out.log",
            err = config['base_log_outdir']+"/{prefix}/cgmlstfinder/run/{sample}_err.log"
        threads:
            8
        conda:
            "../envs/cgmlstfinder.yaml"
        shell:
            """
            mkdir -p {output.cgmlstfinder_out}
            resources/tools/cgmlstfinder/cgMLST.py -s ecoli -db /projects/AusGEM/databases/cgmlstfinder_db -o {output.cgmlstfinder_out} -k $CONDA_PREFIX/bin/kma {input.r1_filt} {input.r2_filt} 
            cp {output.cgmlstfinder_out}/ecoli_summary.txt {output.cgmlst_summary}
            """