import re
import subprocess
import os

if path.exists("resources/tools/pmlst") == False:
    print('pmlst directory not located, downloading pmlst...')
    os.system("git clone https://git@bitbucket.org/genomicepidemiology/pmlst.git resources/tools/pmlst")
    os.system("git clone https://git@bitbucket.org/genomicepidemiology/pmlst_db.git resources/tools/pmlst/pmlst_db")


prefix = config['prefix']
maxthreads = snakemake.utils.available_cpu_count()
scheme = config['pmlst_scheme']


logs = config['base_log_outdir']

rule pMLST_run:
    input:
        assembly = config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta"
    output:
        config['outdir']+"/{prefix}/pMLST/{scheme}/{sample}.out/results.txt"
    conda:
        "../envs/pMLST.yaml"
    params:
        output_dir = config['outdir']+"/{prefix}/pMLST/{scheme}/{sample}.out",
        tmp = config['outdir']+"/{prefix}/pMLST/{scheme}/temp_{sample}",
        pmlst_tool = config['pmlst_tool'],
        pmlst_script_path = config['pMLST_script_path'],
        pmlst_db_path = config['pMLST_db_path'],
        db = "{scheme}",
        scheme = scheme,
        blastn_path = config['blast_bin']

    shell:
        """
        python3 {params.pmlst_script_path} -i {input} -o {params.output_dir} -p {params.pmlst_db_path} {params.pmlst_tool} {params.blastn_path}/blastn -s {params.db} -x -t {params.tmp}
        rm -rf {params.tmp}
        """

rule name_files_pMLST:
    input:
        config['outdir']+"/{prefix}/pMLST/{scheme}/{sample}.out/results.txt"
    output:
        config['outdir']+"/{prefix}/pMLST/{scheme}/{sample}.out/results_named.txt"
    shell:
        """
        awk 'NR == 1 {{print "name\t" $0; next;}}{{print FILENAME "\t" $0;}}' {input} > {output}
        """
