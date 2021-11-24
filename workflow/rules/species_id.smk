import re
import subprocess
import os

prefix = config['prefix']
maxthreads = snakemake.utils.available_cpu_count()

if path.exists("resources/tools/Bracken/bracken") == False:
    print('Pointfinder directory not located, downloading pointfinder...')
    os.system("git clone https://github.com/jenniferlu717/Bracken.git resources/tools/Bracken")

if config['input_type'] == 'assemblies':
    rule run_kraken2:
        input:
            db = config['krakendb'],
            assembly = config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta"
        output:
            out = config['outdir']+"/{prefix}/kraken2/{sample}.out",
            report = config['outdir']+"/{prefix}/kraken2/{sample}.report"
        log:
            config['base_log_outdir']+"/{prefix}/kraken2/run/{sample}_err.log"
        conda:
            "../envs/kraken2.yaml"
        shell:
            """
            kraken2 --db {input.db} --use-names --report {output.report} --output {output.out} {input.assembly}  2> {log}
            """

elif config['input_type'] == 'raw_reads':
    rule run_kraken2:
        input:
            db = config['krakendb'],
            r1_filt = config['outdir']+"/{prefix}/fastp/{sample}.R1.fastq.gz",
            r2_filt = config['outdir']+"/{prefix}/fastp/{sample}.R2.fastq.gz"
        output:
            out = config['outdir']+"/{prefix}/kraken2/{sample}.out",
            report = config['outdir']+"/{prefix}/kraken2/{sample}.report"
        log:
            config['base_log_outdir']+"/{prefix}/kraken2/run/{sample}_err.log"
        conda:
            "../envs/kraken2.yaml"
        shell:
            """
            kraken2 --db {input.db} --use-names --report {output.report} --output {output.out} {input.r1_filt} {input.r2_filt}  2> {log}
            """
            
rule bracken:
    input:
        config['outdir']+"/{prefix}/kraken2/{sample}.report"
    output:
        #Assembly statistics
        #assembly_stats=config['output_dir']+"/assembly-stats/{assembler}/{sample}.txt"
        bracken=config['outdir']+"/{prefix}/bracken/{sample}.bracken.txt",
        species=config['outdir']+"/{prefix}/bracken/{sample}_bracken_species_report.txt"
    params:
        krakendb=config['krakendb']
        #extra="-t",
    log:
        config['base_log_outdir']+"/{prefix}/bracken/{sample}.bracken.log",
    threads: 4
    conda:
        "../envs/kraken2.yaml"
    shell:
        "resources/tools/Bracken/bracken -d {params.krakendb} -i {input} -o {output.bracken} -w {output.species} -r 100 -l S -t {threads} 2>&1 {log}"
    #wrapper:
        #"0.2.0/bio/assembly-stats"
    #    "https://raw.githubusercontent.com/maxlcummins/snakemake-wrappers/assembly-stats/bio/assembly-stats/wrapper.py"

rule run_bracken_summarise:
    input:
        #assembly-stats text files
        bracken_reports=expand(config['outdir']+"/{prefix}/bracken/{sample}.bracken.txt", prefix=prefix, sample=sample_ids)
    output:
        #combined, tab-delimited assembly statistics
        combined_bracken=config['outdir']+"/{prefix}/summaries/bracken_report.txt"
    params:
        extra="",
    log:
        "logs/{prefix}/summaries/combine_bracken.log",
    threads: 1
    script:
        "../../scripts/combine_bracken.py"
