import re
import subprocess
import os

prefix = config['prefix']
maxthreads = snakemake.utils.available_cpu_count()
scheme = config['pmlst_scheme']


logs = config['base_log_outdir']

rule dfast_db_download:
    output:
        directory("resources/dbs/dfast")
    conda:
        "../envs/dfast.yaml"
    shell:
        """
        if [ ! -d "resources/dbs/dfast/protein" ]; then echo 'dfast protein dbs directory not located, downloading to resources/dbs/dfast/protein...'
        #git clone https://github.com/nigyta/dfast_core.git
        dfast_file_downloader.py --dbroot resources/dbs/dfast --protein dfast
        fi
        if [ ! -d "resources/dbs/dfast/hmm" ]; then echo 'dfast hmm dbs directory not located, downloading to resources/dbs/dfast/hmm...'
        dfast_file_downloader.py --dbroot resources/dbs/dfast --cdd Cog --hmm TIGR
        fi
        """

rule dfast_run:
    input:
        assembly = config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta",
        database = "resources/dbs/dfast"
    output:
        directory(config['outdir']+"/{prefix}/dfast/{sample}.out")
    conda:
        "../envs/dfast.yaml"
    params:
        dfast_out = config['outdir']+"/{prefix}/dfast/{sample}.out",
    shell:
        """
        dfast -g {input.assembly} --dbroot {input.database} -o {output} --use_locustag_as_gene_id
        """

rule gff_rename:
    input:
        config['outdir']+"/{prefix}/dfast/{sample}.out"
    output:
        config['outdir']+"/{prefix}/dfast/gffs/{sample}.gff"
    conda:
        "../envs/dfast.yaml"
    shell:
        """
        mv {input}/genome.gff {output}
        """





#rule name_files_pMLST:
#    input:
#        config['outdir']+"/{prefix}/pMLST/{scheme}/{sample}.out/results.txt"
#    output:
#        config['outdir']+"/{prefix}/pMLST/{scheme}/{sample}.out/results_named.txt"
#    shell:
#        """
#        awk 'NR == 1 {{print "name\t" $0; next;}}{{print FILENAME "\t" $0;}}' {input} > {output}
#        """
