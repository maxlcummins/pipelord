import re
import subprocess
import os

if path.exists("resources/dbs/mobsuite/ncbi_plasmid_full_seqs.fas.msh") == False:
    os.system("mkdir -p resources/dbs/mobsuite")
    rule mobsuite_init:
        output:
            "resources/dbs/mobsuite/ncbi_plasmid_full_seqs.fas.msh"
        conda:
            "../envs/mobsuite.yaml"
        log:
            err = config['base_log_outdir']+"/"+config['prefix']+"/mobsuite/init/init.err.log",
            out = config['base_log_outdir']+"/"+config['prefix']+"/mobsuite/init/init.out.log"
        threads:
            1
        shell:
            """
            mob_init -d resources/dbs/mobsuite 1> {log.out} 2> {log.err}
            """

rule mobsuite_recon:
    input:
        assembly = config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta",
        db = "resources/dbs/mobsuite/ncbi_plasmid_full_seqs.fas.msh"
    output:
        config['outdir']+"/{prefix}/mobsuite/{sample}/contig_report.txt"
    log:
        err = config['base_log_outdir']+"/{prefix}/mobsuite/recon/{sample}.err.log",
        out = config['base_log_outdir']+"/{prefix}/mobsuite/recon/{sample}.out.log"
    conda:
        "../envs/mobsuite.yaml"
    threads:
        1
    params:
        outdir = config['outdir']+"/{prefix}/mobsuite/{sample}"
    shell:
        "mob_recon --force --num_threads {threads} --infile {input.assembly} --outdir {params.outdir} 1> {log.out} 2> {log.err}"
        
rule mobsuite_summarise:
    input:
        mobsuite = expand(config['outdir']+"/{prefix}/mobsuite/{sample}/contig_report.txt", prefix=prefix, sample=sample_ids)
    output:
        config['outdir']+"/{prefix}/summaries/mobtyper_results.txt"
    log:
        err = config['base_log_outdir']+"/{prefix}/mobsuite/recon/summarise.err.log",
        out = config['base_log_outdir']+"/{prefix}/mobsuite/recon/summarise.out.log"
    conda:
        "../envs/abritamr.yaml"
    threads:
        1
    script:
        "../../scripts/combine_mobsuite.py"