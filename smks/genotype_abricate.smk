import re
import subprocess
import os

configfile: "misc/masterconfig.yaml"

# Get assemblies
sample_ids, = glob_wildcards(config['raw_reads_path']+"/{sample}.R1.fastq.gz")
prefix = config['prefix']
maxthreads = snakemake.utils.available_cpu_count()

db_location = config['gene_db_location']
gene_dbs = expand(config['gene_dbs'])

logs = config['base_log_outdir']

#print(sample_ids)
#print(db_location)
#print(gene_dbs)


#rule all:
#    input:
#        expand(config['outdir']+"/{prefix}/abricate/{gene_db}/{sample}.tab",
#               sample=sample_ids, gene_db=gene_dbs, prefix=prefix),
#        expand(config['outdir']+"/{prefix}/summaries/abricate_hits.txt", prefix=prefix)
    #expand(config['outdir']+"/{prefix}/summaries/{gene_db}_hits.txt", gene_db=gene_dbs, prefix=prefix)

rule abricate_run:
    input:
        assembly = config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta"
    output:
        tab = config['outdir']+"/{prefix}/abricate/{gene_db}/{sample}.tab",
    log:
        config['base_log_outdir']+"/{prefix}/abricate/run/{gene_db}/{sample}.log"
    conda:
        "config/abricate.yaml"
    threads:
        6
    params:
        db = "{gene_db}",
        gene_db = gene_dbs,
        datadir = config['gene_db_location']
    shell:
        "abricate --nopath --datadir {params.datadir} --db {params.db} {input.assembly} > {output}"

# in beta
rule concatenate_abricate_hits:
    input:
        expand(config['outdir']+"/{prefix}/abricate/{gene_db}/{sample}.tab", sample=sample_ids, gene_db=gene_dbs, prefix=prefix)
    output:
        config['outdir']+"/{prefix}/summaries/abricate_hits.txt"
    log:
        config['base_log_outdir']+"/{prefix}/abricate/concatenate.log"
    conda:
        "config/abricate.yaml"
    threads:
        maxthreads
    shell:
        """
        cat {input} | awk 'BEGIN{{FS=OFS=";"}} {{gsub(/\.fasta/, "", $1)}} 1' | sed -E '1s/#FILE/name/g' | grep -v "#FILE" > {output}
        """

#rule abricate_summary:
#    input:
#        expand(config['outdir']+"/{prefix}/abricate/{gene_db}/{sample}.tab", sample=sample_ids, gene_db=gene_dbs, prefix=prefix)
#    output:
#        config['outdir']+"/{prefix}/summaries/{gene_db}_hits.txt"
#    log:
#        config['base_log_outdir']+"/{prefix}/abricate/{gene_db}_summarise.log"
#    conda:
#        "config/abricate.yaml"
#    threads:
#        maxthreads
#    shell:
#        """
#       abricate --summary {input} > {output}
#        """



#include: "genome_assembly.smk"