import re
import subprocess
import os

# Assumes config is already defined and loaded
prefix = config['prefix']
scheme = config['pmlst_scheme']
logs = config['base_log_outdir']

if config["genotype_modules"]['genome_annotater'] == "dfast":

    rule dfast_db_download:
        output:
            directory("resources/dbs/dfast")
        conda:
            "../envs/dfast.yaml"
        threads: 1
        shell:
            """
            if [ ! -d "resources/dbs/dfast/protein" ]; then 
                echo 'dfast protein dbs directory not located, downloading to resources/dbs/dfast/protein...'
                dfast_file_downloader.py --dbroot resources/dbs/dfast --protein dfast
            fi
            if [ ! -d "resources/dbs/dfast/hmm" ]; then 
                echo 'dfast hmm dbs directory not located, downloading to resources/dbs/dfast/hmm...'
                dfast_file_downloader.py --dbroot resources/dbs/dfast --cdd Cog --hmm TIGR
            fi
            """

    rule dfast_run:
        input:
            assembly = config['outdir'] + "/{prefix}/shovill/assemblies/{sample}.fasta",
            database = "resources/dbs/dfast"
        output:
            directory(config['outdir'] + "/{prefix}/dfast/{sample}.out")
        conda:
            "../envs/dfast.yaml"
        params:
            dfast_out = config['outdir'] + "/{prefix}/dfast/{sample}.out",
        threads: 6
        shell:
            """
            dfast -g {input.assembly} --cpu {threads} --dbroot {input.database} -o {output} --use_locustag_as_gene_id
            """

    rule gff_rename:
        input:
            directory(config['outdir'] + "/{prefix}/dfast/{sample}.out")
        output:
            config['outdir'] + "/{prefix}/dfast/gffs/{sample}.gff"
        conda:
            "../envs/dfast.yaml"
        threads: 1
        shell:
            """
            mv {input}/genome.gff {output}
            """

elif config["genotype_modules"]['genome_annotater'] == "prokka":

    rule prokka_run:
        input:
            assembly = config['outdir'] + "/{prefix}/shovill/assemblies/{sample}.fasta",
        output:
            directory(config['outdir'] + "/{prefix}/prokka/{sample}.out")
        conda:
            "../envs/prokka.yaml"
        params:
            prokka_out = config['outdir'] + "/{prefix}/prokka/{sample}.out",
        threads: 6
        shell:
            """
            prokka --cpus {threads} --outdir {output} {input.assembly} 
            """

    rule prokka_gff_rename:
        input:
            directory(config['outdir'] + "/{prefix}/prokka/{sample}.out")
        output:
            config['outdir'] + "/{prefix}/prokka/gffs/{sample}.gff"
        conda:
            "../envs/prokka.yaml"
        threads: 1
        shell:
            """
            mv {input}/PROKKA_*.gff {output}
            """

elif config["genotype_modules"]['genome_annotater'] == "bakta":

    rule bakta_db_download:
        output:
            directory(config['bakta_db'])
        conda:
            "../envs/bakta.yaml"
        threads: 1
        log:
            config['base_log_outdir'] + "/" + config["prefix"] + "/bakta/bakta_db_dl.log"
        shell:
            """
            if [ ! -d "{output}" ]; then 
                echo 'bakta db not located, downloading to {output}'
                bakta_db download --output {output} > {log} 2>&1
            fi
            """

    rule bakta_run:
        input:
            bakta_db = config['bakta_db'],
            assembly = config['outdir'] + "/{prefix}/shovill/assemblies/{sample}.fasta",
        output:
            directory(config['outdir'] + "/{prefix}/bakta/{sample}.out")
        conda:
            "../envs/bakta.yaml"
        params:
            bakta_db = config['bakta_db'],
            bakta_out = config['outdir'] + "/{prefix}/bakta/{sample}.out",
        threads: 8
        log:
            config['base_log_outdir'] + "/{prefix}/bakta/annotation/{sample}.log"
        shell:
            """
            bakta --db {input.bakta_db} --verbose --output {output} --prefix {wildcards.sample} --threads {threads} {input.assembly}
            # If no .log file is generated, ignore the error with '|| true'
            cp {output}/*.log {log} || true
            """

    rule bakta_gff_rename:
        input:
            directory(config['outdir'] + "/{prefix}/bakta/{sample}.out")
        output:
            config['outdir'] + "/{prefix}/bakta/gffs/{sample}.gff3"
        conda:
            "../envs/bakta.yaml"
        threads: 1
        log:
            config['base_log_outdir'] + "/{prefix}/bakta/rename/{sample}.log"
        shell:
            """
            mv {input}/*.gff3 {output} &> {log}
            """

else:
    print("Please set the genome annotater in the config file to either 'dfast', 'prokka' or 'bakta'. We recommend 'bakta', 'prokka' and 'dfast' in that order of preference.")
