import re
import subprocess
import os
from os import path
import git

#configfile: "misc/masterconfig2.yaml"

# Get assemblies
#sample_ids, = glob_wildcards(config['raw_reads_path']+"/{sample}.R1.fastq.gz")
prefix = config['prefix']
maxthreads = snakemake.utils.available_cpu_count()

if path.exists(config['pointfinder_path']) == False:
    print('Pointfinder directory not located, downloading pointfinder...')
    os.system("git clone https://git@bitbucket.org/genomicepidemiology/pointfinder.git tools/pointfinder")
    os.system("git clone https://git@bitbucket.org/genomicepidemiology/pointfinder_db.git tools/pointfinder/pointfinder_db")



#rule all:
#    input:
#        #expand(config['outdir']+"/{prefix}/pointfinder/{sample}/{sample}_blastn_results_named.tsv", sample=sample_ids, prefix=prefix)
#        expand(config['outdir']+"/{prefix}/summaries/Pointfinder.txt", prefix=prefix)
#        #expand(config['outdir']+"{output_prefix}_abricate.tab", output_prefix = output_pre)

rule pointfinder_run:
    input:
        assembly = config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta"
    output:
        temp(config['outdir']+"/{prefix}/pointfinder/{sample}/{sample}_blastn_results.tsv")
    conda:
        "config/pointfinder.yaml"
    params:
        output_dir = config['outdir']+"/{prefix}/pointfinder/{sample}",
        species = config['pointfinder_species'],
        pointfinder_path =  config['pointfinder_path']
    shell:
        """
        {params.pointfinder_path}/PointFinder.py -i {input} -o {params.output_dir} -p {params.pointfinder_path}/pointfinder_db {params.species} -m blastn -m_p /usr/local/ncbi-blast-ihpc-2.8.1+/bin/blastn
        """

rule name_append:
    input:
        config['outdir']+"/{prefix}/pointfinder/{sample}/{sample}_blastn_results.tsv"
    output:
        temp(config['outdir']+"/{prefix}/pointfinder/{sample}/{sample}_blastn_results_named.tsv")
    shell:
        """awk 'NR == 1 {{print "name\t" $0; next;}}{{print FILENAME "\t" $0;}}' {input} > {output}"""


rule combine:
    input:
        expand(config['outdir']+"/{prefix}/pointfinder/{sample}/{sample}_blastn_results_named.tsv", sample=sample_ids, prefix=prefix)
    output:
        temp(config['outdir']+"/{prefix}/pointfinder/Pointfinder_temp.txt")
    threads:
        maxthreads
    shell:
        """cat {input} > {output}"""

rule clean:
    input:
        config['outdir']+"/{prefix}/pointfinder/Pointfinder_temp.txt"
    output:
        config['outdir']+"/{prefix}/summaries/Pointfinder.txt"
    threads:
        maxthreads
    shell:
        """
        awk 'FNR==1 {{ header = $0; print }} $0 != header' {input} > {output}
        perl -p -i -e 's@.*/@@g' {output}
        perl -p -i -e 's@_blastn_results.tsv@@g' {output}

        """


# rule data_combine:
#	input:
#		config['outdir']+"{sample}/{sample}_blastn_results_named.tsv"
#	output:
#		config['outdir']+"{sample}.dummy_file"
#	params:
#		outfile = config['outdir']+"pointfinder.txt",
#	shell:
#		"""
#		touch {params.outfile}
#		awk 'NR == 1 {{print $0 "\tname_file"; next;}}{{print $0 "\t" FILENAME;}}' {input} >> {params.outdir}
#		touch {output}
#		"""
