import re
import subprocess
import os

prefix = config['prefix']

db_location = config['gene_db_location']
gene_dbs = expand(config['gene_dbs'])
plasmid_screen_db = expand(config['plasmid_screen_db'])

abritamr_outdir = config["outdir"]+"/"+config['prefix']+"/abritamr",

logs = config['base_log_outdir']

rule update_database:
    output:
        temp(config['outdir']+"/{prefix}/abritamr/dummy_out.txt")
    conda:
        "../envs/abritamr.yaml"
    threads:
        1
    shell:
        """
        amrfinder -U
        touch {output}
        """

rule abritamr_single_fofn:
    input:
        config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta"
    output:
        temp(config['outdir']+"/{prefix}/abritamr/fofns/separate/{sample}_fofn.txt")
    conda:
        "../envs/abritamr.yaml"
    threads:
        1
    params:
        out_prefix=config["outdir"]+"/{prefix}/abritamr/"
    shell:
        """
        sample_name=$(basename {input} | sed 's/\.fasta//')
        sample_path=$(realpath --relative-to="." {input})
        echo -e {params.out_prefix}${{sample_name}}'\t'${{sample_path}} >> {output}
        """

rule abritamr_fofn:
    input:
        expand(config['outdir']+"/{prefix}/abritamr/fofns/separate/{sample}_fofn.txt", prefix=prefix, sample=sample_ids)
    output:
        temp(config['outdir']+"/{prefix}/abritamr/fofns/full_fofn.txt")
    conda:
        "../envs/abritamr.yaml"
    threads:
        1
    shell:
        """
        cat {input} > {output}
        """

rule run_abritamr:
    input:
        fofn = config['outdir']+"/{prefix}/abritamr/fofns/full_fofn.txt",
        dummy_db_file = config['outdir']+"/{prefix}/abritamr/dummy_out.txt"
    output:
        config["outdir"]+"/{prefix}/summaries/abritamr_summary_matches.txt",
        config["outdir"]+"/{prefix}/summaries/abritamr_summary_partials.txt",
        config["outdir"]+"/{prefix}/summaries/abritamr_summary_virulence.txt"
    params:
        abritamr_out = abritamr_outdir,
        species = config['abritamr_species']
    log:
        "logs/{prefix}/summaries/combine_abritamr.log"
    threads: 50
    conda:
        "../envs/abritamr.yaml"
    shell:
        """
        abritamr run --prefix {params.abritamr_out} --contigs {input.fofn}  --jobs {threads} --species {params.species}
        mv summary_matches.txt {output[0]}
        mv summary_partials.txt {output[1]}
        mv summary_virulence.txt {output[2]}
        """

ruleorder: update_database > abritamr_single_fofn > abritamr_fofn > run_abritamr