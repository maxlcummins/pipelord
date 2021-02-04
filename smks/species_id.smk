import re
import subprocess
import os

configfile: "misc/masterconfig.yaml"

#print(assembly_path)

# Get reads
sample_ids, = glob_wildcards(config['raw_reads_path']+"/{sample}.R1.fastq.gz")
prefix = config['prefix']
maxthreads = snakemake.utils.available_cpu_count()

#print(config['krakendb'])
#print(sample_ids)

# one rule to rule them all
#rule all:
#    input:
#        expand(config['outdir']+"/{prefix}/kraken2/{sample}.out", sample=sample_ids, prefix=prefix),
#        expand(config['outdir']+"/{prefix}/kraken2/{sample}.report",sample=sample_ids, prefix=prefix),
#        expand(config['outdir']+"/{prefix}/summaries/kraken2_summary.txt", prefix=prefix),


rule run_kraken2:
    input:
        db = config['krakendb'],
        assembly = config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta"
    output:
        out = config['outdir']+"/{prefix}/kraken2/{sample}.out",
        report = config['outdir']+"/{prefix}/kraken2/{sample}.report"
    log:
        config['base_log_outdir']+"/{prefix}/kraken2/run/{sample}.log"
    conda:
        "config/kraken2.yaml"
    threads: maxthreads
    shell:
        """
        kraken2 --db {input.db} --use-names --report {output.report} --output {output.out} {input.assembly}
        """

rule name_kraken2_inline:
    input:
        config['outdir']+"/{prefix}/kraken2/{sample}.report"
    output:
        config['outdir']+"/{prefix}/kraken2/{sample}_clean.report"

    shell:
        """
        awk 'NR == 1 {{print "name_file\t" $0; next;}}{{print FILENAME "\t" $0;}}' {input} > {output}
        perl -p -i -e 's@\.report@@g' {output}
        """

rule combine_kraken2_reports:
    input:
        expand(config['outdir']+"/{prefix}/kraken2/{sample}_clean.report",sample=sample_ids, prefix=prefix)
    output:
        summary_temp = config['outdir']+"/{prefix}/summaries/kraken2_summary_temp.txt",
        summary = config['outdir']+"/{prefix}/summaries/kraken2_summary.txt"

    shell:
        """
        cat {input} > {output.summary_temp}
        #Removes duplicate headers (in this case lines starting with filename)
        awk 'FNR==1 {{ header = $0; print }} $0 != header' {output.summary_temp} > {output.summary}
        """

include: "genome_assembly.smk"