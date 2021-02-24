import re
import subprocess
import os

if path.exists(config['pMLST_db_path']) == False:
    print('pmlst directory not located, downloading pmlst...')
    os.system("git clone https://git@bitbucket.org/genomicepidemiology/pmlst.git tools/pmlst")
    os.system("git clone https://git@bitbucket.org/genomicepidemiology/pmlst_db.git tools/pmlst/pmlst_db")

#include: "genome_assembly.smk"

#configfile: "misc/masterconfig2.yaml"

#Get assemblies
#sample_ids, = glob_wildcards(config['raw_reads_path']+"/{sample}.R1.fastq.gz")
prefix = config['prefix']
maxthreads = snakemake.utils.available_cpu_count()
scheme = config['pmlst_scheme']

#print(sample_ids)

logs = config['base_log_outdir']

#rule all:
#    input:
#        expand(config['outdir']+"/{prefix}/pMLST/{scheme}/{sample}.out/results_tab_named.tsv", sample=sample_ids, prefix=prefix, scheme=scheme),
#        expand(config['outdir']+"/{prefix}/pMLST/{scheme}_pMLST.txt", prefix=prefix, scheme=scheme)

rule pMLST_run:
    input:
        assembly = config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta"
    output:
        config['outdir']+"/{prefix}/pMLST/{scheme}/{sample}.out/results.txt"
    conda:
        "config/pMLST.yaml"
    params:
        output_dir = config['outdir']+"/{prefix}/pMLST/{scheme}/{sample}.out",
        tmp = config['outdir']+"/{prefix}/pMLST/{scheme}/temp_{sample}",
        pmlst_tool = config['pmlst_tool'],
        pmlst_script_path = config['pMLST_script_path'],
        pmlst_db_path = config['pMLST_db_path'],
        db = "{scheme}",
        scheme = scheme,

    shell:
        """
        python3 {params.pmlst_script_path} -i {input} -o {params.output_dir} -p {params.pmlst_db_path} {params.pmlst_tool} -s {params.db} -x -t {params.tmp}
        rm -rf {params.tmp}
        """

# Define input files
#def pmlst_combine_inputs(wildcards):
#    files = expand(config['outdir']+"/{prefix}/pMLST/{scheme}/{sample}.out/results.txt", sample=sample_ids, prefix=prefix, scheme=scheme)
#    return files

rule name_files_pMLST:
    input:
        config['outdir']+"/{prefix}/pMLST/{scheme}/{sample}.out/results.txt"
    output:
        config['outdir']+"/{prefix}/pMLST/{scheme}/{sample}.out/results_named.txt"
    shell:
        """
        awk 'NR == 1 {{print "name\t" $0; next;}}{{print FILENAME "\t" $0;}}' {input} > {output}
        """

#rule list_pmlst_outs:
#    input:
#        expand(config['outdir']+"/{prefix}/pMLST/{scheme}/{sample}.out/results_named.txt", sample=sample_ids, prefix=prefix, scheme=scheme)
#        #['a,b,c']
#    output:
#        config['outdir']+"/{prefix}/summaries/pMLST_tmp.txt"
#    params:
#        outdir=config['outdir']+"/{prefix}/pMLST"
#    shell:
#        """
#        scripts/test.sh -i {input} -o {output}
#        """
#
#rule pmlst_clean:
#    threads:
#        maxthreads
#    input:
#        config['outdir']+"/{prefix}/summaries/pMLST_tmp.txt"
#    output:
#        config['outdir']+"/{prefix}/summaries/pMLST.txt"
#    #log:
#        #config['base_log_outdir']+"/{prefix}/pMLST/{scheme}/combine/combine.log"
#    threads:
#        1
#    shell:
#        """
#        grep -E "pMLST profile|Sequence Type" {input} | perl -p -e 's@pMLST profile: @@g' | awk '{{printf "%s%s",$0,NR%2?"\t":RS}}' > {output}
#        #Trim the first column to generate a column with pMLST scheme
#        perl -p -i -e 's@^.*Inc@@g' {output}
#        #Trim the first column to generate a column with pMLST scheme
#        perl -p -i -e 's@^.*pBSSB1@pBSSB1@g' {output}
#        #Trim column two to generate a column with sample names
#        perl -p -i -e 's@\t.*(inc[^/]+|pbssb1[^/]+)/@\t@g' {output}
#        # Trim column with sample names to clean them
#        perl -p -i -e 's@.out/results.txt@@g' {output}
#        #Trim pMLST column to get rid of junk text
#        perl -p -i -e 's@Sequence Type: @@g' {output}
#        #Clean square brackets from pMLST names
#        perl -p -i -e 's@(\]|\[)@@g' {output}
#        #Fix scheme column names
#        perl -p -i -e 's@^@Inc@g' {output}
#        perl -p -i -e 's@Inc @Inc@g' {output}
#        perl -p -i -e 's@^.*family@pbssb1-family@g'  {output}
#        #Change order of columns
#        awk -F'\t' -v OFS="\t" '{{ print $2, $1, $3}}' {output} > tmp && mv tmp {output}
#        #Insert a header
#        echo -e "name\tpMLST_scheme\tpMLST" | cat - {output} > tmp && mv tmp {output}
#        #Remove MLST scheme/sample combinations with no hits
#        grep -v 'Unknown' {output} > tmp && mv tmp {output}
#        """