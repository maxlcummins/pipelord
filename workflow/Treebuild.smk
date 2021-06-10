import re
import subprocess
import os
from os import path
import git

configfile: "config/citrobacter_all.yaml"

#Number of threads
maxthreads = snakemake.utils.available_cpu_count()

outdir = config['outdir']
prefix = config['prefix']
gene_dbs = expand(config['gene_dbs'])
plasmid_screen_db = expand(config['plasmid_screen_db'])

scheme = config['pmlst_scheme']
email = config['email']


subset_fofn = config['subset_fofn']
subset_names = config['subset_name']
with open(subset_fofn) as f:
        subset = [line.strip("\n") for line in f]


rule all:
    input:
        expand(config['outdir']+"/{prefix}/roary/{subset_name}", prefix=prefix, subset_name=subset_names),
        expand(config['outdir']+"/{prefix}/snp_sites/{subset_name}/core_gene_alignment_snp_sites.aln", prefix=prefix, subset_name=subset_names),
        expand(config['outdir']+"/{prefix}/trees/{subset_name}/CGA_snp_sites.treefile", prefix=prefix, subset_name=subset_names),
        expand(config['outdir']+"/{prefix}/trees/{subset_name}/CGA_full.treefile", prefix=prefix, subset_name=subset_names),
        expand(config['outdir']+"/{prefix}/snp_dists/{subset_name}/CGA_snp_sites_pairwise_snps_counts.csv", prefix=prefix, subset_name=subset_names),
        expand(config['outdir']+"/{prefix}/snp_dists/{subset_name}/CGA_full_pairwise_snps_counts.csv", prefix=prefix, subset_name=subset_names),


rule roary_subset:
    input:
        config['outdir']+"/{prefix}/dfast/gffs/{subset}.gff"
    output:
        config['outdir']+"/{prefix}/subsets/{subset_name}/gffs/{subset}.gff"
    shell:
        """
        cp -n {input} {output}
        """

rule roary:
    input:
        expand(config['outdir']+"/{prefix}/subsets/{subset_name}/gffs/{subset}.gff", prefix=prefix, subset_name=subset_names, subset=subset)
    output:
        directory(config['outdir']+"/{prefix}/roary/{subset_name}")
    conda:
        "envs/roary.yaml"
    log:
        config['base_log_outdir']+"/{prefix}/roary/{subset_name}.log"
    threads:
        maxthreads
    params:
        out = config['outdir']+"/{prefix}/roary/{subset_name}",
        core_req = config['roary_core_req'],
        kraken_db = config['krakendb']
    shell:
        "roary -e -v -n -r -cd {params.core_req} -qc -k {params.kraken_db} -f {params.out} {input}  2> {log}"

rule snp_sites:
    input:
        config['outdir']+"/{prefix}/roary/{subset_name}"
    output:
        config['outdir']+"/{prefix}/snp_sites/{subset_name}/core_gene_alignment_snp_sites.aln"
    conda:
        "envs/snp_sites.yaml"
    log:
        config['base_log_outdir']+"/{prefix}/roary/logs/{subset_name}/snp_sites.log"
    shell:
        "snp-sites -c {input}/core_gene_alignment.aln > {output}"

rule iqtree:
    input:
        config['outdir']+"/{prefix}/roary/{subset_name}"
    output:
        config['outdir']+"/{prefix}/trees/{subset_name}/CGA_full.treefile"
    conda:
        "envs/iqtree.yaml"
    log:
        config['base_log_outdir']+"/{prefix}/roary/logs/{subset_name}/iqtree.log"
    params:
        out = config['outdir']+"/{prefix}/trees/{subset_name}/CGA_full"
    shell:
        """
        iqtree -s {input}/core_gene_alignment.aln -pre {params.out} -m MFP -bb 1000
        """

rule iqtree_snp_sites:
    input:
        config['outdir']+"/{prefix}/snp_sites/{subset_name}/core_gene_alignment_snp_sites.aln"
    output:
        config['outdir']+"/{prefix}/trees/{subset_name}/CGA_snp_sites.treefile"
    conda:
        "envs/iqtree.yaml"
    log:
        config['base_log_outdir']+"/{prefix}/roary/logs/{subset_name}/iqtree_snp_sites.log"
    params:
        out = config['outdir']+"/{prefix}/trees/{subset_name}/CGA_snp_sites"
    shell:
        """
        iqtree -s {input} -pre {params.out} -m MFP -bb 1000
        """

rule snp_dists_snp_sites:
        input:
                config['outdir']+"/{prefix}/snp_sites/{subset_name}/core_gene_alignment_snp_sites.aln"
        output:
                config['outdir']+"/{prefix}/snp_dists/{subset_name}/CGA_snp_sites_pairwise_snps_counts.csv"
        conda:
                "envs/snp_dists.yaml"
        shell:
                "snp-dists -c {input} > {output}"



rule snp_dists_full:
        input:
                config['outdir']+"/{prefix}/roary/{subset_name}"
        output:
                config['outdir']+"/{prefix}/snp_dists/{subset_name}/CGA_full_pairwise_snps_counts.csv"
        conda:
                "envs/snp_dists.yaml"
        shell:
                "snp-dists -c {input}/core_gene_alignment.aln > {output}"


onerror:
    print("An error occurred")
    #Email user
    if config['email'] != '':
        shell("echo Your Snakemake job with prefix \'{prefix}\' has failed. | mail -s 'Snakemake job has failed' {email}")


# Please code deities forgive my frankenstein onsuccess script. I will repackage it in a future distribution into a script... for now I just need it to run!
onsuccess:
#Announce job finish
    print("Job finished")
    #Email user
    if config['email'] != '':
        shell("echo Your Snakemake job with prefix \'{prefix}\' has finished. It has been written to \'{outdir}/{prefix}/summaries/\' | mail -s 'Snakemake job has finished' {email}")


ruleorder: roary > snp_sites > iqtree > iqtree_snp_sites > snp_dists_snp_sit