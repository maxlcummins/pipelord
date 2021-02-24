import re
import subprocess
import os

prefix = config['prefix']
maxthreads = snakemake.utils.available_cpu_count()


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

rule name_kraken2_inline:
    input:
        config['outdir']+"/{prefix}/kraken2/{sample}.report"
    output:
        config['outdir']+"/{prefix}/kraken2/{sample}_clean.report"
    shell:
        """
        awk 'NR == 1 {{print FILENAME "\t" $0; next;}}{{print FILENAME "\t" $0;}}' {input} > {output}
        perl -p -i -e 's@\.report@@g' {output}
        """
