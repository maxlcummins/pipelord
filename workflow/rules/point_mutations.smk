import re
import subprocess
import os
from os import path
import git

prefix = config['prefix']
maxthreads = snakemake.utils.available_cpu_count()

if path.exists("resources/tools/pointfinder") == False:
    print('Pointfinder directory not located, downloading pointfinder...')
    os.system("git clone https://git@bitbucket.org/genomicepidemiology/pointfinder.git resources/tools/pointfinder")
    os.system("git clone https://git@bitbucket.org/genomicepidemiology/pointfinder_db.git resources/tools/pointfinder/pointfinder_db")
    os.system("perl -p -i -e 's@^sample_name = filename.split.*@sample_name = filename.split(\".\")[0]@g' resources/tools/pointfinder/PointFinder.py")

rule pointfinder_run:
    input:
        assembly = config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta"
    output:
        config['outdir']+"/{prefix}/pointfinder/{sample}/{sample}_blastn_results.tsv"
    conda:
        "../envs/pointfinder.yaml"
    log:
        config['base_log_outdir']+"/{prefix}/pointfinder/pointfinder_run/{sample}.log"
    params:
        output_dir = config['outdir']+"/{prefix}/pointfinder/{sample}",
        species = config['pointfinder_species'],
        blastn_path = config['blast_bin']
    shell:
        """
        resources/tools/pointfinder/PointFinder.py -i {input} -o {params.output_dir} -p resources/tools/pointfinder/pointfinder_db {params.species} -m blastn -m_p {params.blastn_path}/blastn 2> {log}
        """

rule name_append:
    input:
        config['outdir']+"/{prefix}/pointfinder/{sample}/{sample}_blastn_results.tsv"
    output:
        config['outdir']+"/{prefix}/pointfinder/{sample}/{sample}_blastn_results_named.tsv"
    shell:
        """awk 'NR == 1 {{print "name\t" $0; next;}}{{print FILENAME "\t" $0;}}' {input} > {output}"""

